///module with small utility fucntions

//use bitvec::prelude::*;
//use packed_seq::{PackedSeq, Seq};

///hashing function for u64 using xorshift
pub fn xorshift_u64(mut x: u64) -> u64 {
    x ^= x<<13;
    x ^= x>>7;
    x ^= x<<17;
    x
}

pub fn _xorshift_u32(mut x: u32) -> u32 {
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    x
}

///since no implementation directly on u128 exists, I just used whatever numbers of shifts
pub fn xorshift_u128(mut x: u128) -> u128 {
    x ^= x << 17;
    x ^= x >> 23;
    x ^= x << 5;
    x
}

///returns the amount of positive bits in a boolean vector
pub fn sum_vec_bool(boolean_vector: &Vec<bool>) -> usize {
    let mut counter: usize = 0;
    for val in boolean_vector {
        if *val {counter += 1};
    }
    counter
}

#[cfg(test)]
mod tests {
    use super::{_xorshift_u32, sum_vec_bool, xorshift_u64, xorshift_u128};

    #[test]
    fn xorshift_u64_is_deterministic() {
        assert_eq!(xorshift_u64(42), xorshift_u64(42));
        assert_ne!(xorshift_u64(42), 42);
    }

    #[test]
    fn xorshift_u128_is_deterministic() {
        assert_eq!(xorshift_u128(42), xorshift_u128(42));
        assert_ne!(xorshift_u128(42), 42);
    }

    #[test]
    fn sum_vec_bool_counts_true_values() {
        let values = vec![true, false, true, true, false];
        assert_eq!(sum_vec_bool(&values), 3);
    }

    #[test]
    fn xorshift_u64_zero_stays_zero() {
        assert_eq!(xorshift_u64(0), 0);
    }

    #[test]
    fn xorshift_u128_zero_stays_zero() {
        assert_eq!(xorshift_u128(0), 0);
    }

    #[test]
    fn xorshift_u32_is_deterministic() {
        assert_eq!(_xorshift_u32(12345), _xorshift_u32(12345));
        assert_ne!(_xorshift_u32(12345), 12345);
    }

    #[test]
    fn sum_vec_bool_empty_vector_is_zero() {
        let values = vec![];
        assert_eq!(sum_vec_bool(&values), 0);
    }

    #[test]
    fn sum_vec_bool_all_true_counts_every_entry() {
        let values = vec![true, true, true, true];
        assert_eq!(sum_vec_bool(&values), 4);
    }
}
