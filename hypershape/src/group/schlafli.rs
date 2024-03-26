use hypermath::collections::{approx_hashmap::ApproxHashMapKey, ApproxHashMap};
use hypermath::prelude::*;
use itertools::Itertools;

use super::{CoxeterGroup, GroupError, GroupResult, IsometryGroup};

/// Schlafli symbol for a convex polytope.
#[derive(Debug, Default, Clone, PartialEq, Eq, Hash)]
pub struct SchlafliSymbol {
    indices: Vec<Vec<usize>>,
}
impl SchlafliSymbol {
    /// Constructs an integer Schlafli symbol.
    pub fn from_linear_indices(indices: Vec<usize>) -> Self {
        let dim = indices.len() + 1;
        Self {
            indices: (0..dim)
                .map(|i| {
                    (0..dim)
                        .map(|j| {
                            if i == j {
                                1
                            } else if i == j + 1 || j == i + 1 {
                                indices[std::cmp::min(i, j)]
                            } else {
                                2
                            }
                        })
                        .collect_vec()
                })
                .collect_vec(),
        }
    }

    /// Constructs and validates integer Schlafli symbol.
    pub fn from_matrix_indices(indices: Vec<Vec<usize>>) -> GroupResult<Self> {
        let dim = indices.len();
        // Index matrix is not square
        if indices.iter().any(|r| r.len() != dim) {
            return Err(GroupError::CDInvalid);
        }

        for i in 0..dim {
            for j in 0..=i {
                // Index matrix has non-positive integer
                if indices[i][j] < 1 {
                    return Err(GroupError::CDInvalid);
                }
                // Index matrix has ones off diagonal or no ones on diagonal
                if (i == j) != (indices[i][j] == 1) {
                    return Err(GroupError::CDInvalid);
                }
                // Index matrix is not symmetric
                if indices[i][j] != indices[j][i] {
                    return Err(GroupError::CDInvalid);
                }
            }
        }
        Ok(Self { indices })
    }

    /// Convert Coxeter group to SchlafliSymbol.
    pub fn from_coxeter_group(c: CoxeterGroup) -> SchlafliSymbol {
        Self {
            indices: (0..c.generator_count())
                .map(|i| {
                    (0..c.generator_count())
                        .map(|j| c.coxeter_matrix_element(i, j) as usize)
                        .collect_vec()
                })
                .collect_vec(),
        }
    }

    /// Returns the indices of the Schlafli symbol.
    pub fn indices(&self) -> &[Vec<usize>] {
        &self.indices
    }

    /// Constructs an integer Schlafli symbol from a string.
    pub fn from_string(string: &str) -> Self {
        let xs = string
            .split(',')
            .map(|s| s.trim().parse().unwrap_or(0))
            .collect_vec();
        Self::from_linear_indices(xs)
    }

    /// Number of dimensions of the polytope described by the Schlafli symbol.
    pub fn ndim(&self) -> u8 {
        self.indices.len() as u8
    }

    /// Returns the i,j-th entry of the Schläfli matrix.
    pub fn mirror_dot(&self, i: usize, j: usize) -> f64 {
        -(std::f64::consts::PI as Float / self.indices[i][j] as Float).cos()
    }

    /// Returns the list of mirrors.
    pub fn mirrors(&self) -> GroupResult<Vec<Mirror>> {
        let mut ret = vec![];
        // The final mirror vectors will look like this, with each row as a
        // vector:
        //
        // ```
        // [ ? 0 0 0 0 ]
        // [ ? ? 0 0 0 ]
        // [ ? ? ? 0 0 ]
        // [ ? ? ? ? 0 ]
        // [ ? ? ? ? ? ]
        // ```
        //
        // If this matrix is `L`, `L Lᵀ = A`, where `A` is the Schläfli
        // matrix of the Coxeter-Dynkin diagram. This is a Cholesky
        // decomposition. We use the Cholesky–Banachiewicz algorithm.
        // https://en.wikipedia.org/wiki/Cholesky_decomposition#Computation
        for i in 0..self.indices.len() {
            ret.push(Mirror(Vector::zero(self.ndim())));
            for j in 0..=i {
                let mut sum = 0.0;
                for k in 0..j {
                    sum += ret[i].0[k as u8] * ret[j].0[k as u8];
                }

                if i == j {
                    let val = self.mirror_dot(i, i) - sum;
                    if val < 0.0 {
                        return Err(GroupError::CDHyperbolic);
                    }
                    ret[i].0[j as u8] = val.sqrt();
                } else {
                    ret[i].0[j as u8] = 1.0 / ret[j].0[j as u8] * (self.mirror_dot(i, j) - sum);
                }
            }
        }
        Ok(ret)
    }

    /// Returns a matrix that transforms from the mirror basis (where each
    /// component of the vector gives a distance from a mirror plane) to the
    /// base space.
    pub fn mirror_basis(&self) -> GroupResult<Matrix> {
        Matrix::from_cols(self.mirrors()?.into_iter().map(|Mirror(m)| m))
            .transpose()
            .inverse()
            .ok_or(GroupError::CDEuclidean)
    }

    /// Returns the list of mirrors as generators.
    pub fn generators(&self) -> GroupResult<Vec<Isometry>> {
        Ok(self.mirrors()?.into_iter().map(|m| m.into()).collect())
    }

    /// Constructs the isometry group described by the Schlafli symbol.
    pub fn group(&self) -> GroupResult<IsometryGroup> {
        IsometryGroup::from_generators(&self.generators()?)
    }

    /// Expands an object by the symmetry.
    pub fn expand<T: Clone + ApproxHashMapKey>(
        &self,
        object: T,
        transform: fn(&Isometry, &T) -> T,
    ) -> GroupResult<Vec<T>> {
        let generators = self.generators()?;

        let mut seen = ApproxHashMap::new();
        seen.insert(&object, ());

        let mut next_unprocessed_index = 0;
        let mut ret = vec![object];
        while next_unprocessed_index < ret.len() {
            let unprocessed_object = ret[next_unprocessed_index].clone();
            for gen in &generators {
                let new_object = transform(gen, &unprocessed_object);
                if seen.insert(&new_object, ()).is_none() {
                    ret.push(new_object);
                }
            }
            next_unprocessed_index += 1;
        }
        Ok(ret)
    }
}

#[derive(Debug, Clone, PartialEq)]
struct MirrorGenerator {
    mirrors: Vec<Mirror>,
}
impl From<MirrorGenerator> for Isometry {
    fn from(gen: MirrorGenerator) -> Self {
        gen.mirrors
            .into_iter()
            .map(Isometry::from)
            .fold(Isometry::ident(), |a, b| a * b)
    }
}
impl From<MirrorGenerator> for Matrix {
    fn from(gen: MirrorGenerator) -> Self {
        gen.mirrors
            .into_iter()
            .map(Matrix::from)
            .fold(Matrix::EMPTY_IDENT, |a, b| a * b)
    }
}

/// Mirror hyperplane that intersects the origin, defined by its normal vector.
#[derive(Debug, Clone, PartialEq)]
pub struct Mirror(pub Vector);
impl From<Mirror> for Isometry {
    fn from(Mirror(v): Mirror) -> Self {
        Isometry::from_reflection_normalized(v)
    }
}
impl From<Mirror> for Matrix {
    fn from(Mirror(v): Mirror) -> Self {
        let ndim = v.ndim();
        let mut ret = Matrix::ident(ndim);
        for x in 0..ndim {
            for y in 0..ndim {
                *ret.get_mut(x, y) = ret.get(x, y) - 2.0 * v[x] * v[y];
            }
        }
        ret
    }
}

#[cfg(test)]
mod tests {
    use super::super::Group;
    use super::*;

    #[test]
    fn test_cube_group() {
        let g = SchlafliSymbol::from_linear_indices(vec![4, 3])
            .group()
            .unwrap();

        assert_eq!(48, g.element_count());
    }
}
