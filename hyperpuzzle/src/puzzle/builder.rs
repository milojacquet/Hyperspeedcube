#![allow(clippy::too_many_arguments, clippy::too_many_lines)]

use std::collections::{hash_map, HashMap, HashSet};
use std::ops::Range;
use std::sync::{Arc, Weak};

use eyre::{bail, ensure, OptionExt, Result, WrapErr};
use hypermath::prelude::*;
use hypershape::prelude::*;
use itertools::Itertools;
use parking_lot::Mutex;
use smallvec::smallvec;

use super::simplices::{Simplexifier, VertexId};
use super::{
    centroid::Centroid, Color, ColorInfo, Facet, FacetSet, Mesh, MeshVertexData, Notation,
    PerColor, PerFacet, PerPiece, PerSticker, Piece, PieceInfo, PieceSet, PieceType, PieceTypeInfo,
    Puzzle, StickerInfo,
};

/// Puzzle being constructed.
#[derive(Debug)]
pub struct PuzzleBuilder {
    /// Puzzle ID.
    pub id: String,
    /// Name of the puzzle.
    pub name: String,

    /// Space where the puzzle exists.
    pub space: Arc<Mutex<Space>>,
    /// Puzzle pieces.
    pub pieces: PerPiece<PieceBuilder>,
    /// Sticker colors.
    pub colors: PerColor<ColorBuilder>,
}
impl PuzzleBuilder {
    /// Returns a builder for a puzzle that starts with an empty Euclidean space
    /// with no pieces.
    pub fn new_soup(name: String, id: String, ndim: u8) -> Result<Self> {
        Ok(PuzzleBuilder {
            name,
            id,
            space: Arc::new(Mutex::new(Space::new(ndim)?)),
            pieces: PerPiece::new(),
            colors: PerColor::new(),
        })
    }
    /// Returns a builder for a puzzle that starts with a single solid piece
    /// occupying all of Euclidean space.
    pub fn new_solid(name: String, id: String, ndim: u8) -> Result<(Self, Piece)> {
        let mut this = PuzzleBuilder::new_soup(name, id, ndim)?;
        this.pieces.push(PieceBuilder {
            shape: this.space.lock().whole_space(),
            is_active: true,
        })?;
        Ok((this, Piece(0)))
    }

    /// Cuts each piece by a cut, throwing away the portions that are outside
    /// the cut. Every piece in the old set becomes inactive, and each piece in
    /// the new set inherits its active status from the corresponding piece in
    /// the old set.
    pub fn carve(&mut self, pieces: &PieceSet, cut_manifold: ManifoldRef) -> Result<PieceSet> {
        let mut cut = AtomicCut::carve(cut_manifold);
        self.cut_and_deactivate_pieces(&mut cut, pieces)
    }

    /// Cuts each piece by a cut, keeping all results. Every piece in the old
    /// set becomes inactive, and each piece in the new set inherits its active
    /// status from the corresponding piece in the old set.
    pub fn slice(&mut self, pieces: &PieceSet, cut_manifold: ManifoldRef) -> Result<PieceSet> {
        let mut cut = AtomicCut::slice(cut_manifold);
        self.cut_and_deactivate_pieces(&mut cut, pieces)
    }

    /// Returns the set of active pieces.
    pub fn active_pieces(&self) -> PieceSet {
        self.pieces
            .iter()
            .filter(|(_id, piece)| piece.is_active)
            .map(|(id, _piece)| id)
            .collect()
    }

    fn cut_and_deactivate_pieces(
        &mut self,
        cut: &mut AtomicCut,
        pieces: &PieceSet,
    ) -> Result<PieceSet> {
        let mut space = self.space.lock();

        let old_pieces = pieces;
        let mut new_pieces = PieceSet::new();
        for piece in old_pieces.iter() {
            let old_piece = &mut self.pieces[piece];

            // Cut and deactivate piece.
            for new_piece in cut_and_deactivate_piece(&mut space, old_piece, cut)? {
                let new_piece_id = self.pieces.push(new_piece)?;
                new_pieces.insert(new_piece_id);
            }
        }

        Ok(new_pieces)
    }

    /// Adds a new color, assigned to a manifold.
    pub fn add_color(&mut self, manifold: ManifoldRef) -> Result<Color> {
        Ok(self.colors.push(ColorBuilder {
            manifold,
            default_color: None,
            name: None,
        })?)
    }
    /// Sets the name for a color.
    pub fn set_color_name(&mut self, color: Color, name: String) -> Result<()> {
        self.colors
            .get_mut(color)
            .ok_or_eyre("index out of range")?
            .name = Some(name);
        Ok(())
    }
    /// Sets the default for a color.
    pub fn set_color_default_color(&mut self, color: Color, default_color: String) -> Result<()> {
        self.colors
            .get_mut(color)
            .ok_or_eyre("index out of range")?
            .default_color = Some(default_color);
        Ok(())
    }
    /// Sets the order of all colors, given a list of the new color order. Each
    /// `i`th element in `new_order` is the ID of the old color that should be
    /// the new `i`th color.
    pub fn set_color_order(&mut self, new_order: PerColor<Color>) -> Result<()> {
        let mut old_colors = std::mem::take(&mut self.colors).map(|_, color| Some(color));
        self.colors = new_order.try_map(|_, old_color_id| {
            old_colors
                .get_mut(old_color_id)
                .ok_or_eyre("index out of range")?
                .take()
                .ok_or_eyre("duplicate color in order")
        })?;
        Ok(())
    }

    /// Returns a map from manifold to color.
    fn manifold_colors(&self) -> HashMap<ManifoldRef, Color> {
        self.colors
            .iter()
            .map(|(color_id, color_builder)| (color_builder.manifold, color_id))
            .collect()
    }
    /// Returns a list of `ColorInfo`s to be used in the final `Puzzle`.
    fn color_infos(&self) -> PerColor<ColorInfo> {
        let used_names: HashSet<String> = self
            .colors
            .iter_values()
            .filter_map(|color| color.default_color.clone())
            .collect();
        let mut unused_names = uppercase_names().filter(|name| !used_names.contains(name));
        self.colors.map_ref(|_, color| ColorInfo {
            name: (color.name.clone())
                .unwrap_or_else(|| unused_names.next().expect("ran out of names")),
            default_color: color.default_color.clone(),
        })
    }

    /// Performs the final steps of building a puzzle, generating the mesh and
    /// assigning IDs to pieces, stickers, etc.
    pub fn build(self) -> Result<Arc<Puzzle>> {
        let space = self.space.lock();

        let mut mesh = Mesh::new_empty(space.ndim());
        mesh.color_count = self.colors.len();

        // Only colored manifolds have an entry in `manifold_colors`.
        let manifold_colors: HashMap<ManifoldRef, Color> = self.manifold_colors();
        // All manifolds have an entry in `manifold_to_facet`.
        let mut manifold_to_facet: HashMap<ManifoldRef, Facet> = HashMap::new();

        // As we construct the mesh, we'll renumber all the pieces and stickers
        // to exclude inactive ones.
        let mut pieces = PerPiece::<PieceInfo>::new();
        let mut stickers = PerSticker::<StickerInfo>::new();
        let mut facets = PerFacet::<TempFacetData>::new();
        let colors = self.color_infos();

        let mut simplexifier = Simplexifier::new(&space);

        // Construct the mesh for each piece.
        for (_old_piece_id, piece) in &self.pieces {
            // Skip inactive pieces.
            if !piece.is_active {
                continue;
            }

            let piece_centroid = simplexifier.shape_centroid_point(piece.shape.id)?;

            let piece_id = pieces.push(PieceInfo {
                stickers: smallvec![],
                piece_type: PieceType(0), // TODO: piece types
                centroid: piece_centroid.clone(),
            })?;

            // Add stickers to mesh, sorted ordered by color. It's important
            // that internal stickers are processed last, so that they are all
            // in a consecutive range for `piece_internals_index_ranges`.
            let mut piece_stickers = space
                .boundary_of(piece.shape)
                .map(|sticker_shape| {
                    let facet_manifold = space.manifold_of(sticker_shape);
                    let color = *manifold_colors
                        .get(&facet_manifold)
                        .unwrap_or(&Color::INTERNAL);
                    (color, sticker_shape)
                })
                // Skip internals for 4D+.
                .filter(|(color, _sticker_shape)| space.ndim() < 4 || *color != Color::INTERNAL)
                .collect_vec();
            piece_stickers.sort();

            let sticker_shrink_vectors = compute_sticker_shrink_vectors(
                &space,
                &mut simplexifier,
                piece.shape,
                &piece_stickers,
            )?;

            let mut piece_internals_indices_start = None;

            for (sticker_color, sticker_shape) in piece_stickers {
                if sticker_color != Color::INTERNAL {
                    let sticker_id = stickers.push(StickerInfo {
                        piece: piece_id,
                        color: sticker_color,
                    })?;
                    pieces[piece_id].stickers.push(sticker_id);
                }

                let sticker_centroid = simplexifier.shape_centroid(sticker_shape.id)?;
                let manifold = space.manifold_of(sticker_shape);
                let facet_id = match manifold_to_facet.entry(manifold) {
                    hash_map::Entry::Occupied(e) => *e.get(),
                    hash_map::Entry::Vacant(e) => {
                        let facet_id = facets.push(TempFacetData::new(&space, manifold)?)?;
                        *e.insert(facet_id)
                    }
                };

                facets[facet_id].centroid += sticker_centroid;

                let (polygon_index_range, triangles_index_range, edges_index_range) =
                    build_shape_polygons(
                        &space,
                        &mut mesh,
                        &mut simplexifier,
                        &sticker_shrink_vectors,
                        sticker_shape,
                        piece_id,
                        facet_id,
                    )?;

                if sticker_color == Color::INTERNAL {
                    if piece_internals_indices_start.is_none() {
                        piece_internals_indices_start = Some((
                            polygon_index_range.start,
                            triangles_index_range.start,
                            edges_index_range.start,
                        ));
                    }
                } else {
                    mesh.add_sticker(
                        polygon_index_range,
                        triangles_index_range,
                        edges_index_range,
                    )?;
                }
            }

            let mut piece_internals_polygon_range = 0..0;
            let mut piece_internals_triangle_range = 0..0;
            let mut piece_internals_edge_range = 0..0;
            if let Some((polygon_start, tri_start, edge_start)) = piece_internals_indices_start {
                piece_internals_polygon_range = polygon_start..mesh.polygon_count;
                piece_internals_triangle_range = tri_start..mesh.triangle_count() as u32;
                piece_internals_edge_range = edge_start..mesh.edge_count() as u32;
            }
            mesh.add_piece(
                &piece_centroid,
                piece_internals_polygon_range,
                piece_internals_triangle_range,
                piece_internals_edge_range,
            )?;
        }

        for (_, facet_data) in facets {
            mesh.add_facet(facet_data.centroid.center(), facet_data.normal)?;
        }

        Ok(Arc::new_cyclic(|this| Puzzle {
            this: Weak::clone(this),
            name: self.name,
            id: self.id,

            mesh,

            pieces,
            stickers,
            piece_types: [PieceTypeInfo {
                name: "Piece".to_string(), // TODO piece types
            }]
            .into_iter()
            .collect(),
            colors,

            scramble_moves_count: 500, // TODO

            notation: Notation {},
        }))
    }
}

/// Piece of a puzzle during puzzle construction.
#[derive(Debug, Clone)]
pub struct PieceBuilder {
    /// Shape of the piece.
    pub shape: AtomicPolytopeRef,
    /// Whether the piece should be part of the final puzzle.
    pub is_active: bool,
}

/// Sticker color during puzzle construction.
#[derive(Debug, Clone)]
pub struct ColorBuilder {
    /// Manifold of the color; stickers flush with this manifold will be
    /// assigned this color.
    pub manifold: ManifoldRef,

    /// Name for the color, which will be automatically chosen if omitted.
    pub name: Option<String>,
    /// Default color string.
    pub default_color: Option<String>,
}

#[derive(Debug)]
struct TempFacetData {
    centroid: Centroid,
    normal: Vector,
}
impl TempFacetData {
    fn new(space: &Space, manifold: ManifoldRef) -> Result<Self> {
        let ipns_blade = space.blade_of(manifold).opns_to_ipns(space.ndim());
        ensure!(
            ipns_blade.ipns_is_flat(),
            "4D backface culling assumes flat faces",
        );

        Ok(Self {
            centroid: Centroid::ZERO,
            normal: ipns_blade
                .ipns_plane_normal()
                .ok_or_eyre("no plane normal")?,
        })
    }
}

/// Computes the sticker shrink vectors for a piece.
///
/// Each vertex shrinks along a vector pointing toward the centroid of the
/// piece, projected onto whatever sticker facets the vertex is part of. For
/// example, if a vertex is on an edge (1D manifold) of a 3D polytope, then its
/// shrink vector will point toward the centroid of the piece, projected onto
/// that edge. If a vertex is on a corner of its polytope, then its shrink
/// vector is zero.
fn compute_sticker_shrink_vectors(
    space: &Space,
    simplexifier: &mut Simplexifier<'_>,
    piece_shape: AtomicPolytopeRef,
    stickers: &[(Color, AtomicPolytopeRef)],
) -> Result<HashMap<VertexId, Vector>> {
    // For the purposes of sticker shrink, we don't care about internal facets.
    let colored_sticker_shapes = stickers
        .iter()
        .filter(|&&(color, _sticker_shape)| color != Color::INTERNAL)
        .map(|&(_color, sticker_shape)| sticker_shape)
        .collect_vec();

    let ndim = space.ndim();

    // Make our own facet IDs that will stay within this object. We only care
    // about the facets that this piece has stickers on.
    let mut facet_blades_ipns = PerFacet::new();

    #[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
    enum Element {
        Point(VertexId),
        NonPoint(AtomicPolytopeId),
    }
    fn sum_centroids(
        simplexifier: &mut Simplexifier<'_>,
        elements: impl Iterator<Item = Element>,
    ) -> Result<Option<Centroid>> {
        let mut ret = Centroid::ZERO;
        for element in elements {
            ret += match element {
                Element::Point(p) => Ok(Centroid::new(&simplexifier[p], 1.0)),
                Element::NonPoint(p) => simplexifier.shape_centroid(p),
            }?;
        }
        Ok((!ret.is_zero()).then_some(ret))
    }

    // For each element of the polytope, compute a set of the facet manifolds
    // that have a sticker containing the element.
    let mut facet_set_per_vertex: HashMap<VertexId, FacetSet> = HashMap::new();
    let mut elements_and_facet_sets_by_rank: Vec<HashMap<Element, FacetSet>> =
        vec![HashMap::new(); space.ndim() as usize + 1];
    for &sticker_shape in &colored_sticker_shapes {
        let manifold = space.manifold_of(sticker_shape);
        let facet = facet_blades_ipns.push(space.blade_of(manifold).opns_to_ipns(ndim))?;

        for vertex in simplexifier.vertex_set(sticker_shape)? {
            facet_set_per_vertex
                .entry(vertex)
                .or_default()
                .insert(facet);
        }
        for element in space.elements_of(sticker_shape.id)? {
            if space.ndim_of(element) > 0 {
                elements_and_facet_sets_by_rank[space.ndim_of(element) as usize]
                    .entry(Element::NonPoint(element))
                    .or_default()
                    .insert(facet);
            }
        }
    }
    elements_and_facet_sets_by_rank[0] = facet_set_per_vertex
        .iter()
        .map(|(&vertex, facet_set)| (Element::Point(vertex), facet_set.clone()))
        .collect();

    // Find the largest (by rank) elements contained by all the sticker facets
    // of the piece.
    let centroid_of_greatest_common_elements: Option<Centroid> = elements_and_facet_sets_by_rank
        .iter()
        .rev()
        .map(|elements_and_facet_sets| {
            // Find elements that are contained by all sticker facets of the
            // piece.
            let elements_with_maximal_facet_set = elements_and_facet_sets
                .iter()
                .filter(|(_element, facet_set)| facet_set.len() == colored_sticker_shapes.len())
                .map(|(element, _facet_set)| *element);
            // Add up their centroids. Technically we should take the centroid
            // of their convex hull, but this works well enough.
            sum_centroids(simplexifier, elements_with_maximal_facet_set)
        })
        // Select the elements with the largest rank.
        .find_map(|result_option| result_option.transpose())
        .transpose()?;
    // If such elements exist, then all vertices can shrink to the same point.
    if let Some(centroid) = centroid_of_greatest_common_elements {
        let shrink_target = centroid.center();
        return Ok(simplexifier
            .vertex_set(piece_shape)?
            .into_iter()
            .map(|vertex| {
                let vertex_position = &simplexifier[vertex];
                (vertex, &shrink_target - vertex_position)
            })
            .collect());
    }

    // Otherwise, find the best elements for each set of facets. If a vertex is
    // not contained by any facets, then it will shrink toward the centroid of
    // the piece.
    let piece_centroid = simplexifier.shape_centroid_point(piece_shape.id)?;

    // Compute the shrink target for each possible facet set that has a good
    // shrink target.
    let unique_facet_sets_of_vertices = elements_and_facet_sets_by_rank[0].values().unique();
    let shrink_target_by_facet_set: HashMap<&FacetSet, Vector> = unique_facet_sets_of_vertices
        .map(|facet_set| {
            // Find the largest elements of the piece that are contained by all
            // the facets in this set. There must be at least one vertex.
            let centroid_of_greatest_common_elements: Centroid = elements_and_facet_sets_by_rank
                .iter()
                .rev()
                .map(|elements_and_facet_sets| {
                    // Find elements that are contained by all sticker facets of
                    // the vertex.
                    let elements_with_superset_of_facets = elements_and_facet_sets
                        .iter()
                        .filter(|(_element, fs)| facet_set.iter().all(|f| fs.contains(f)))
                        .map(|(element, _fs)| *element);
                    // Add up their centroids. Technically we should take the
                    // centroid of their convex hull, but this works well
                    // enough.
                    sum_centroids(simplexifier, elements_with_superset_of_facets)
                })
                // Select the elements with the largest rank.
                .find_map(|result_option| result_option.transpose())
                // There must be some element with a superset of `facet_set`
                // because `facet_set` came from a vertex.
                .expect("no element with facet subset")?;

            eyre::Ok((facet_set, centroid_of_greatest_common_elements.center()))
        })
        .try_collect()?;

    // Compute shrink vectors for all vertices.
    let shrink_vectors = simplexifier
        .vertex_set(piece_shape)?
        .into_iter()
        .map(|vertex| {
            let facet_set = facet_set_per_vertex.remove(&vertex).unwrap_or_default();
            let vertex_pos = &simplexifier[vertex];
            let shrink_vector = match shrink_target_by_facet_set.get(&facet_set) {
                Some(target) => target - vertex_pos,
                None => &piece_centroid - vertex_pos,
            };

            (vertex, shrink_vector)
        });
    Ok(shrink_vectors.collect())
}

fn cut_and_deactivate_piece(
    space: &mut Space,
    piece: &mut PieceBuilder,
    cut: &mut AtomicCut,
) -> Result<Vec<PieceBuilder>> {
    // Deactivate old piece.
    let is_active = std::mem::replace(&mut piece.is_active, false);

    Ok(space
        .cut_atomic_polytope_set([piece.shape].into_iter().collect(), cut)
        .context("error cutting piece")?
        .into_iter()
        .map(|shape| PieceBuilder { shape, is_active })
        .collect())
}

fn build_shape_polygons(
    space: &Space,
    mesh: &mut Mesh,
    simplexifier: &mut Simplexifier<'_>,
    sticker_shrink_vectors: &HashMap<VertexId, Vector>,
    sticker_shape: AtomicPolytopeRef,
    piece_id: Piece,
    facet_id: Facet,
) -> Result<(Range<usize>, Range<u32>, Range<u32>)> {
    let polygons_start = mesh.polygon_count;
    let triangles_start = mesh.triangle_count() as u32;
    let edges_start = mesh.edge_count() as u32;

    for polygon in space.children_with_ndim(sticker_shape, 2) {
        let polygon_id = mesh.next_polygon_id()?;

        // Get the tangent space so that we can compute tangent vectors
        // for each vertex.
        let manifold = space.manifold_of(polygon);
        let blade = space.blade_of(manifold);
        ensure!(blade.grade() == 4, "polygon must lie on 2D manifold");
        let tangent_space = blade.opns_tangent_space();

        // Triangulate the polygon.
        let tris = simplexifier.triangles(polygon)?;

        // The simplexifier and mesh each have their own set of vertex IDs, so
        // we need to be able to map between them.
        let mut vertex_id_map: HashMap<VertexId, u32> = HashMap::new();
        for old_vertex_ids in tris {
            let mut new_vertex_ids = [0; 3];
            for (i, old_vertex_id) in old_vertex_ids.into_iter().enumerate() {
                new_vertex_ids[i] = match vertex_id_map.entry(old_vertex_id) {
                    hash_map::Entry::Occupied(e) => *e.get(),
                    hash_map::Entry::Vacant(e) => {
                        let position = &simplexifier[old_vertex_id];
                        let tangent_vectors = tangent_space.at(position);
                        let Some([u_tangent, v_tangent]) = tangent_vectors.as_deref() else {
                            bail!("bad tangent space");
                        };

                        let sticker_shrink_vector = sticker_shrink_vectors
                            .get(&old_vertex_id)
                            .ok_or_eyre("missing sticker shrink vector for vertex")?;

                        let new_vertex_id = mesh.add_vertex(MeshVertexData {
                            position,
                            u_tangent,
                            v_tangent,
                            sticker_shrink_vector,
                            piece_id,
                            facet_id,
                            polygon_id,
                        });
                        *e.insert(new_vertex_id)
                    }
                };
            }
            mesh.triangles.push(new_vertex_ids);
        }

        for edge in simplexifier.polygon_edges(polygon)? {
            // We should have seen all these vertices before because they show
            // up in triangles.
            mesh.edges.push(edge.map(|id| vertex_id_map[&id]))
        }
    }

    let polygons_end = mesh.polygon_count;
    let triangles_end = mesh.triangle_count() as u32;
    let edges_end = mesh.edge_count() as u32;
    Ok((
        polygons_start..polygons_end,
        triangles_start..triangles_end,
        edges_start..edges_end,
    ))
}

fn uppercase_names() -> impl Iterator<Item = String> {
    fn string_from_base_26(bytes: &[u8]) -> String {
        bytes.iter().rev().map(|&byte| byte as char).collect()
    }

    let mut digits = vec![];
    std::iter::from_fn(move || {
        for char in &mut digits {
            *char += 1;
            if *char <= b'Z' {
                return Some(string_from_base_26(&digits));
            }
            *char = b'A';
        }
        digits.push(b'A');
        Some(string_from_base_26(&digits))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_uppercase_names() {
        let mut iter = uppercase_names();
        assert_eq!(iter.next().as_deref(), Some("A"));
        assert_eq!(iter.next().as_deref(), Some("B"));
        assert_eq!(iter.next().as_deref(), Some("C"));
        let mut iter = iter.skip(22);
        assert_eq!(iter.next().as_deref(), Some("Z"));
        assert_eq!(iter.next().as_deref(), Some("AA"));
        assert_eq!(iter.next().as_deref(), Some("AB"));
        let mut iter = iter.skip(26);
        assert_eq!(iter.next().as_deref(), Some("BC"));
        let mut iter = iter.skip(645);
        assert_eq!(iter.next().as_deref(), Some("ZY"));
        assert_eq!(iter.next().as_deref(), Some("ZZ"));
        assert_eq!(iter.next().as_deref(), Some("AAA"));
        assert_eq!(iter.next().as_deref(), Some("AAB"));
    }
}
