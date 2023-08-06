//! Common mathematical utility functions that didn't fit anywhere else.

use std::ops::{Add, Mul};

use super::Float;

/// Linearly interpolates (unclamped) between two values.
pub fn mix<T>(a: T, b: T, t: Float) -> <T::Output as Add>::Output
where
    T: Mul<Float>,
    T::Output: Add,
{
    a * (1.0 - t) + b * t
}

/// Returns the element of an iterator with the minimum value, allowing floats
/// or other `PartialOrd` types.
pub fn min_by_key<T, K: PartialOrd>(
    elems: impl IntoIterator<Item = T>,
    mut f: impl FnMut(&T) -> K,
) -> Option<T> {
    let mut iter = elems.into_iter();
    let mut min_elem = iter.next()?;
    let mut min_key = f(&min_elem);
    for elem in iter {
        let key = f(&elem);
        if key < min_key {
            min_elem = elem;
            min_key = key;
        }
    }
    Some(min_elem)
}

/// Divides `lhs` by `rhs` if the reciprocal of `rhs` is finite; otherwise
/// returns `None`.
pub fn try_div<T>(lhs: T, rhs: Float) -> Option<T::Output>
where
    T: Mul<Float>,
{
    let recip_rhs = rhs.recip();
    recip_rhs.is_finite().then(|| lhs * recip_rhs)
}

/// Returns the square root of `n` if the result is finite; otherwise returns
/// `None`.
pub fn try_sqrt(n: Float) -> Option<Float> {
    let ret = n.sqrt();
    ret.is_finite().then_some(ret)
}

/// Iterator with a manually-specified exact size.
#[derive(Debug, Clone)]
pub struct WithExactSizeIter<I> {
    iter: I,
    len: usize,
}
impl<I: Iterator> Iterator for WithExactSizeIter<I> {
    type Item = I::Item;

    fn next(&mut self) -> Option<Self::Item> {
        if self.len > 0 {
            self.len -= 1;
        }
        self.iter.next()
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.len, Some(self.len))
    }
}
impl<I: Iterator> ExactSizeIterator for WithExactSizeIter<I> {
    fn len(&self) -> usize {
        self.len
    }
}

/// Extension trait for `.with_exact_size()`.
pub trait IterWithExactSizeExt: Iterator + Sized {
    /// Returns an `ExactSizeIterator` that thinks it has length `len`.
    ///
    /// This length is not checked.
    fn with_exact_size(self, len: usize) -> WithExactSizeIter<Self>;
}
impl<I: Iterator> IterWithExactSizeExt for I {
    fn with_exact_size(self, len: usize) -> WithExactSizeIter<Self> {
        WithExactSizeIter { iter: self, len }
    }
}

/// Stolen from
/// https://github.com/rust-lang/rust/blob/e6ce5627a9e8af9ae4673a390954fffaf526e5cc/library/core/src/num/int_macros.rs#L2204-L2222
///
/// When #![feature(int_roundings)] is merged, delete this.
pub fn next_multiple_of(lhs: u64, rhs: u64) -> u64 {
    let m = lhs % rhs;

    if m == 0 {
        lhs
    } else {
        lhs + (rhs - m)
    }
}
