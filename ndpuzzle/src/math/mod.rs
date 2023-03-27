//! N-dimensional vector math library.

pub use approx::AbsDiffEq;
use num_traits::Zero;

#[macro_use]
mod impl_macros;
#[macro_use]
mod vector;
pub mod cga;
mod matrix;
pub mod permutations;
mod sign;
pub mod util;

pub use matrix::*;
pub use sign::*;
pub use vector::*;

/// Small floating-point value used for comparisons and tiny offsets.
pub const EPSILON: f32 = 0.0001;

/// Compares two numbers, but considers them equal if they are separated by less
/// than `EPSILON`.
pub fn approx_eq<T: AbsDiffEq<Epsilon = f32>>(a: &T, b: &T) -> bool {
    approx::abs_diff_eq!(a, b, epsilon = EPSILON)
}

/// Compares two numbers, but considers them equal if they are separated by less
/// than `EPSILON`.
pub fn approx_cmp<T: AbsDiffEq<Epsilon = f32> + PartialOrd>(a: &T, b: &T) -> std::cmp::Ordering {
    if approx_eq(a, b) {
        std::cmp::Ordering::Equal
    } else if a < b {
        std::cmp::Ordering::Less
    } else {
        std::cmp::Ordering::Greater
    }
}
/// Returns whether one number is less than another by at least `EPSILON`.
pub fn approx_lt<T: AbsDiffEq<Epsilon = f32> + PartialOrd>(a: &T, b: &T) -> bool {
    a < b && !approx_eq(a, b)
}
/// Returns whether one number is greater than another by at least `EPSILON`.
pub fn approx_gt<T: AbsDiffEq<Epsilon = f32> + PartialOrd>(a: &T, b: &T) -> bool {
    a > b && !approx_eq(a, b)
}

/// Returns whether `x` has an absolute value greater than `EPSILON`.
pub fn is_approx_nonzero<T: AbsDiffEq<Epsilon = f32> + Zero>(x: &T) -> bool {
    !approx_eq(x, &T::zero())
}
/// Returns whether `x` is less than `-EPSILON`.
pub fn is_approx_negative<T: AbsDiffEq<Epsilon = f32> + PartialOrd + Zero>(x: &T) -> bool {
    approx_lt(x, &T::zero())
}
/// Returns whether `x` is greater than `EPSILON`.
pub fn is_approx_positive<T: AbsDiffEq<Epsilon = f32> + PartialOrd + Zero>(x: &T) -> bool {
    approx_gt(x, &T::zero())
}