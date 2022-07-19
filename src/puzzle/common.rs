use itertools::Itertools;
use serde::{Deserialize, Serialize};
use smallvec::SmallVec;
use std::fmt;
use std::ops::{BitOr, BitXorAssign, Index, Mul, Neg};
use strum::{Display, EnumIter, EnumMessage};

use super::*;

#[delegatable_trait]
#[enum_dispatch]
pub trait PuzzleType {
    fn ty(&self) -> PuzzleTypeEnum;
    fn name(&self) -> &str;
    fn family_display_name(&self) -> &'static str;
    fn family_internal_name(&self) -> &'static str;

    fn layer_count(&self) -> u8;
    fn family_max_layer_count(&self) -> u8;

    /// Returns the maximum radius of the puzzle's 3D projection.
    fn projection_radius_3d(&self, p: StickerGeometryParams) -> f32;
    fn scramble_moves_count(&self) -> usize;

    fn faces(&self) -> &[FaceInfo];
    fn pieces(&self) -> &[PieceInfo];
    fn stickers(&self) -> &[StickerInfo];
    fn twist_axes(&self) -> &[TwistAxisInfo];
    fn twist_directions(&self) -> &[TwistDirectionInfo];

    fn twist_axis_from_name(&self, name: &str) -> Option<TwistAxis> {
        (0..self.twist_axes().len() as u8)
            .map(TwistAxis)
            .find(|&twist_axis| self.info(twist_axis).name == name)
    }
    fn twist_direction_from_name(&self, name: &str) -> Option<TwistDirection> {
        (0..self.twist_directions().len() as u8)
            .map(TwistDirection)
            .find(|&twist_direction| self.info(twist_direction).name == name)
    }

    fn check_layers(&self, layers: LayerMask) -> Result<(), &'static str> {
        let layer_count = self.layer_count() as u32;
        if layers.0 > 0 || layers.0 < 1 << layer_count {
            Ok(())
        } else {
            Err("invalid layer mask")
        }
    }
    fn all_layers(&self) -> LayerMask {
        let layer_count = self.layer_count() as u32;
        LayerMask((1 << layer_count) - 1)
    }
    fn slice_layers(&self) -> LayerMask {
        LayerMask((self.all_layers().0 >> 1) & !1)
    }
    fn reverse_layers(&self, layers: LayerMask) -> LayerMask {
        LayerMask(layers.0.reverse_bits() >> (32 - self.layer_count()))
    }

    fn make_recenter_twist(&self, axis: TwistAxis) -> Result<Twist, String>;

    fn reverse_twist(&self, twist: Twist) -> Twist {
        Twist {
            axis: twist.axis,
            direction: self.reverse_twist_direction(twist.direction),
            layers: twist.layers,
        }
    }
    fn canonicalize_twist(&self, twist: Twist) -> Twist;
    fn can_twists_combine(&self, prev: Option<Twist>, curr: Twist, metric: TwistMetric) -> bool;

    fn reverse_twist_direction(&self, direction: TwistDirection) -> TwistDirection;
    fn chain_twist_directions(&self, dirs: &[TwistDirection]) -> Option<TwistDirection>;

    fn twist_command_short_description(
        &self,
        axis_name: Option<TwistAxis>,
        direction: TwistDirection,
        layers: LayerMask,
    ) -> String {
        match axis_name {
            Some(axis) => self.twist_short_description(Twist {
                axis,
                direction,
                layers,
            }),
            None => {
                let dir = self.info(direction).symbol;
                if layers.is_default() {
                    format!("Ø{}", dir)
                } else if layers.is_contiguous_from_outermost() {
                    format!("{}Ø{}", layers.count(), dir)
                } else {
                    format!("{{{}}}Ø{}", layers.long_description(), dir)
                }
            }
        }
    }
    fn twist_short_description(&self, twist: Twist) -> String;
}

trait PuzzleTypeRefExt {
    fn deref_internal(&self) -> Self;
}
#[delegate_to_methods]
#[delegate(PuzzleType, target_ref = "deref_internal")]
impl<'a, P: PuzzleType> PuzzleTypeRefExt for &'a P {
    fn deref_internal(&self) -> &'a P {
        *self
    }
}

#[enum_dispatch]
pub trait PuzzleState: PuzzleType {
    fn twist(&mut self, twist: Twist) -> Result<(), &'static str>;
    fn is_piece_affected_by_twist(&self, twist: Twist, piece: Piece) -> bool {
        twist.layers[self.layer_from_twist_axis(twist.axis, piece)]
    }
    fn pieces_affected_by_twist(&self, twist: Twist) -> Vec<Piece> {
        (0..self.pieces().len() as _)
            .map(Piece)
            .filter(|&piece| self.is_piece_affected_by_twist(twist, piece))
            .collect()
    }
    fn layer_from_twist_axis(&self, twist_axis: TwistAxis, piece: Piece) -> u8;

    fn sticker_geometry(
        &self,
        sticker: Sticker,
        p: StickerGeometryParams,
    ) -> Option<StickerGeometry>;

    fn is_solved(&self) -> bool;
}

/// Enumeration of all puzzle types.
#[derive(Serialize, Deserialize, Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum PuzzleTypeEnum {
    /// 3D Rubik's cube.
    Rubiks3D { layer_count: u8 },
    /// 4D Rubik's cube.
    Rubiks4D { layer_count: u8 },
}
#[delegate_to_methods]
#[delegate(PuzzleType, target_ref = "as_dyn_type")]
impl PuzzleTypeEnum {
    fn as_dyn_type(&self) -> &dyn PuzzleType {
        match *self {
            PuzzleTypeEnum::Rubiks3D { layer_count } => rubiks_3d::puzzle_type(layer_count),
            PuzzleTypeEnum::Rubiks4D { layer_count } => rubiks_4d::puzzle_type(layer_count),
        }
    }
}
impl Default for PuzzleTypeEnum {
    fn default() -> Self {
        Self::Rubiks4D { layer_count: 3 }
    }
}
impl fmt::Display for PuzzleTypeEnum {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.name())
    }
}
impl AsRef<str> for PuzzleTypeEnum {
    fn as_ref(&self) -> &str {
        self.name()
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub struct Twist {
    pub axis: TwistAxis,
    pub direction: TwistDirection,
    pub layers: LayerMask,
}

/// Puzzle of any type.
#[enum_dispatch(PuzzleType, PuzzleState)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Puzzle {
    /// 3D Rubik's cube.
    Rubiks3D(Rubiks3D),
    /// 4D Rubik's cube.
    Rubiks4D(Rubiks4D),
}
impl Default for Puzzle {
    fn default() -> Self {
        Self::new(PuzzleTypeEnum::default())
    }
}
impl Puzzle {
    /// Creates a new puzzle of a particular type.
    pub fn new(ty: PuzzleTypeEnum) -> Puzzle {
        match ty {
            PuzzleTypeEnum::Rubiks3D { layer_count } => {
                Puzzle::Rubiks3D(Rubiks3D::new(layer_count))
            }
            PuzzleTypeEnum::Rubiks4D { layer_count } => {
                Puzzle::Rubiks4D(Rubiks4D::new(layer_count))
            }
        }
    }
}

#[derive(Debug, Default, Copy, Clone, PartialEq, Eq, Hash)]
pub struct Piece(pub u16);
#[derive(Debug, Default, Copy, Clone, PartialEq, Eq, Hash)]
pub struct Sticker(pub u16);
#[derive(Debug, Default, Copy, Clone, PartialEq, Eq, Hash)]
pub struct Face(pub u8);
#[derive(Debug, Default, Copy, Clone, PartialEq, Eq, Hash)]
pub struct TwistAxis(pub u8);
#[derive(Debug, Default, Copy, Clone, PartialEq, Eq, Hash)]
pub struct TwistDirection(pub u8);

pub trait PuzzleInfo<T> {
    type Output;

    fn info(&self, thing: T) -> &Self::Output;
}
macro_rules! impl_puzzle_info_trait {
    (fn $method:ident($thing:ty) -> &$thing_info:ty) => {
        impl<T: PuzzleType + ?Sized> PuzzleInfo<$thing> for T {
            type Output = $thing_info;

            fn info(&self, thing: $thing) -> &$thing_info {
                &self.$method()[thing.0 as usize]
            }
        }
    };
}
impl_puzzle_info_trait!(fn faces(Face) -> &FaceInfo);
impl_puzzle_info_trait!(fn pieces(Piece) -> &PieceInfo);
impl_puzzle_info_trait!(fn stickers(Sticker) -> &StickerInfo);
impl_puzzle_info_trait!(fn twist_axes(TwistAxis) -> &TwistAxisInfo);
impl_puzzle_info_trait!(fn twist_directions(TwistDirection) -> &TwistDirectionInfo);

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct PieceInfo {
    pub stickers: SmallVec<[Sticker; 8]>,
}
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub struct StickerInfo {
    pub piece: Piece,
    pub color: Face,
}
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub struct FaceInfo {
    pub symbol: &'static str, // e.g., "R"
    pub name: &'static str,   // e.g., "Right"
}
impl FaceInfo {
    pub const fn new(symbol: &'static str, name: &'static str) -> Self {
        Self { symbol, name }
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub struct TwistAxisInfo {
    pub name: &'static str, // e.g., "R"
}
impl AsRef<str> for TwistAxisInfo {
    fn as_ref(&self) -> &str {
        self.name
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub struct TwistDirectionInfo {
    pub symbol: &'static str, // "'"
    pub name: &'static str,   // "CCW"
}
impl AsRef<str> for TwistDirectionInfo {
    fn as_ref(&self) -> &str {
        self.name
    }
}
impl TwistDirectionInfo {
    pub const fn new(symbol: &'static str, name: &'static str) -> Self {
        Self { symbol, name }
    }
}

// TODO: revamp turn metrics (see https://www.speedsolving.com/wiki/index.php/Metric)

/// Convention for counting moves.
#[derive(
    Serialize, Deserialize, Debug, Display, EnumIter, EnumMessage, Copy, Clone, PartialEq, Eq, Hash,
)]
#[serde(rename_all = "UPPERCASE")]
pub enum TwistMetric {
    /// Quarter Slice Turn Metric: Each twist counts separately. Whole-puzzle
    /// rotations are not counted.
    #[strum(
        serialize = "QSTM",
        message = "Quarter Slice Turn Metric (default)",
        detailed_message = "Each twist counts separately. Whole-puzzle rotations are not counted."
    )]
    Qstm,

    /// Face Turn Metric: Consecutive twists with the same face and layers are
    /// combined.
    #[strum(
        serialize = "FTM",
        message = "Face Turn Metric",
        detailed_message = "Consecutive twists with the same face and layers are combined."
    )]
    Ftm,

    /// Slice Turn Metric: Consecutive twists with the same face are combined,
    /// even with different layers.
    #[strum(
        serialize = "STM",
        message = "Slice Turn Metric",
        detailed_message = "Consecutive twists with the same face are combined, even with different layers."
    )]
    Stm,

    /// Execution Turn Metric: Each twist counts separately, including
    /// whole-puzzle rotations.
    #[strum(
        serialize = "ETM",
        message = "Execution Turn Metric",
        detailed_message = "Each twist counts separately, including whole-puzzle rotations."
    )]
    Etm,
}
impl Default for TwistMetric {
    fn default() -> Self {
        Self::Qstm
    }
}
impl TwistMetric {
    /// Returns the next twist metric in a cycle.
    pub fn next(self) -> Self {
        match self {
            Self::Qstm => Self::Ftm,
            Self::Ftm => Self::Stm,
            Self::Stm => Self::Etm,
            Self::Etm => Self::Qstm,
        }
    }
}

/// Positive or negative.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Sign {
    /// Positive.
    Pos = 1,
    /// Negative.
    Neg = -1,
}
impl Neg for Sign {
    type Output = Sign;
    fn neg(self) -> Sign {
        match self {
            Sign::Pos => Sign::Neg,
            Sign::Neg => Sign::Pos,
        }
    }
}
impl Mul<Sign> for Sign {
    type Output = Sign;
    fn mul(self, rhs: Sign) -> Sign {
        match self {
            Sign::Pos => rhs,
            Sign::Neg => -rhs,
        }
    }
}
impl Sign {
    /// Returns an integer representation of the sign (either -1 or 1).
    pub const fn int(self) -> isize {
        match self {
            Sign::Pos => 1,
            Sign::Neg => -1,
        }
    }
    /// Returns a floating-point representation of the sign (either -1.0 or
    /// 1.0).
    pub const fn float(self) -> f32 {
        self.int() as f32
    }
    /// Returns an iterator over all signs.
    pub fn iter() -> impl Clone + Iterator<Item = Sign> {
        [Sign::Pos, Sign::Neg].into_iter()
    }
}

/// Bitmask selecting a subset of a puzzle's layers.
#[derive(Serialize, Deserialize, Debug, Copy, Clone, PartialEq, Eq, Hash)]
#[serde(transparent)]
pub struct LayerMask(pub u32);
impl Default for LayerMask {
    fn default() -> Self {
        Self(1)
    }
}
impl Index<u8> for LayerMask {
    type Output = bool;

    fn index(&self, index: u8) -> &Self::Output {
        match self.0 & (1 << index) {
            0 => &false,
            _ => &true,
        }
    }
}
impl LayerMask {
    pub(crate) fn is_default(self) -> bool {
        self == Self::default()
    }
    pub(crate) fn digits(self) -> String {
        // Just give up if there's more than 9 layers.
        (0..9)
            .filter(|&i| self[i])
            .map(|i| (i as u8 + b'1') as char)
            .collect()
    }
    pub(crate) fn short_description(self) -> String {
        match self.count() {
            0 => "none".to_owned(),
            _ => (0..32).filter(|&i| self[i]).map(|i| i + 1).join(", "),
        }
    }
    pub(crate) fn long_description(self) -> String {
        match self.count() {
            0 => "no layers".to_owned(),
            1 => format!("layer {}", self.0.trailing_zeros() + 1),
            _ => format!(
                "layers {}",
                (0..32).filter(|&i| self[i]).map(|i| i + 1).join(", ")
            ),
        }
    }
    pub(crate) fn count(self) -> u32 {
        self.0.count_ones()
    }
    pub(crate) fn is_contiguous_from_outermost(self) -> bool {
        self.0.count_ones() == self.0.trailing_ones()
    }
    pub(crate) fn get_single_layer(self) -> Option<u32> {
        if self.count() == 1 {
            Some(self.0.trailing_zeros())
        } else {
            None
        }
    }
}

/// Bitmasks selecting a subset of a puzzle's twist axes and layers.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub struct TwistSelection {
    /// Bitmask of twist axes.
    pub axis_mask: u32,
    /// Bitmask of layers.
    pub layer_mask: u32,
}
impl Default for TwistSelection {
    fn default() -> Self {
        Self {
            axis_mask: 0_u32,
            layer_mask: 1_u32,
        }
    }
}
impl BitOr for TwistSelection {
    type Output = Self;

    fn bitor(self, rhs: Self) -> Self::Output {
        Self {
            axis_mask: self.axis_mask | rhs.axis_mask,
            layer_mask: self.layer_mask | rhs.layer_mask,
        }
    }
}
impl BitXorAssign for TwistSelection {
    fn bitxor_assign(&mut self, rhs: Self) {
        self.axis_mask ^= rhs.axis_mask;
        self.layer_mask ^= rhs.layer_mask;
    }
}
impl TwistSelection {
    /// Returns the layer mask if any layers are selected, or the default layer
    /// mask (outermost layer only) otherwise.
    pub fn layer_mask_or_default(self) -> LayerMask {
        if self.layer_mask != 0_u32 {
            LayerMask(self.layer_mask)
        } else {
            LayerMask::default()
        }
    }

    /// Returns whether the twist selection includes a particular sticker.
    pub fn has_sticker(self, puzzle: &dyn PuzzleState, sticker: Sticker) -> bool {
        let piece = puzzle.info(sticker).piece;

        // Filter by twist_axis and layer.
        let layer_mask = self.layer_mask_or_default();
        (0..puzzle.twist_axes().len() as _)
            .filter(|i| (self.axis_mask >> i) & 1 != 0)
            .map(TwistAxis)
            .map(|twist_axis| puzzle.layer_from_twist_axis(twist_axis, piece))
            .all(|layer| layer_mask[layer])
    }
}
