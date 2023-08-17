use super::*;

/// ID for an oriented object in a [`Space`].
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub struct SignedRef<I> {
    /// Unoriented ID.
    pub id: I,
    /// Orientation.
    pub sign: Sign,
}
impl<I: fmt::Display> fmt::Display for SignedRef<I> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}", self.sign, self.id)
    }
}
impl<I> From<I> for SignedRef<I> {
    fn from(id: I) -> Self {
        SignedRef {
            id,
            sign: Sign::Pos,
        }
    }
}
impl<I: Fits64> Fits64 for SignedRef<I> {
    unsafe fn from_u64(x: u64) -> Self {
        Self {
            id: I::from_u64(x >> 1),
            sign: if x & 1 == 0 { Sign::Pos } else { Sign::Neg },
        }
    }

    fn to_u64(self) -> u64 {
        (self.id.to_u64() << 1) | self.sign as u64
    }
}
impl<I> Neg for SignedRef<I> {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        self.sign = -self.sign;
        self
    }
}
impl<I: Eq + Fits64> PartialOrd for SignedRef<I> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl<I: Eq + Fits64> Ord for SignedRef<I> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.to_u64().cmp(&other.to_u64())
    }
}
hypermath::impl_mul_sign!(impl<I> Mul<Sign> for SignedRef<I>);
hypermath::impl_mulassign_sign!(impl<I: Clone> MulAssign<Sign> for SignedRef<I>);