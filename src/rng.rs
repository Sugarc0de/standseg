//! glibc-compatible `random(3)`.
//!
//! The C program's tie-break is `#define flip() (random() & 01)`, and it never
//! calls `srandom`, so it runs from glibc's default seed of 1. The golden
//! fixtures were produced on Linux, so reproducing them byte-for-byte means
//! reproducing glibc's generator exactly -- Apple's libc `random()` is a
//! different one, and calling the platform's would silently desync everything
//! downstream. See PLAN.md section 3.1.
//!
//! glibc's default is TYPE_3: a 31-word additive-feedback generator with
//! separation 3, seeded through a Lehmer recurrence and warmed by discarding
//! 10 * 31 outputs.

const DEG: usize = 31;
const SEP: usize = 3;

pub struct GlibcRandom {
    r: [u32; DEG],
    f: usize,
    rear: usize,
}

impl GlibcRandom {
    /// glibc's unseeded state, equivalent to `srandom(1)`.
    pub fn new() -> Self {
        Self::with_seed(1)
    }

    pub fn with_seed(seed: u32) -> Self {
        // glibc maps seed 0 to 1; the state must never be all-zero.
        let seed = if seed == 0 { 1 } else { seed };

        let mut r = [0u32; DEG];
        r[0] = seed;

        // r[i] = (16807 * r[i-1]) % 2147483647, in Schrage form to stay inside
        // 32 bits. glibc does this in signed arithmetic and folds negatives.
        for i in 1..DEG {
            let prev = r[i - 1] as i64;
            let hi = prev / 127773;
            let lo = prev % 127773;
            let mut word = 16807 * lo - 2836 * hi;
            if word < 0 {
                word += 2147483647;
            }
            r[i] = word as u32;
        }

        let mut this = Self {
            r,
            f: SEP,
            rear: 0,
        };

        // Warm up: glibc discards 10 * DEG outputs.
        for _ in 0..(10 * DEG) {
            this.next_u32();
        }
        this
    }

    /// One `random()` result: 31 bits.
    pub fn next_u32(&mut self) -> u32 {
        // Additive feedback with wraparound, then drop the low bit.
        let val = self.r[self.f].wrapping_add(self.r[self.rear]);
        self.r[self.f] = val;
        let result = val >> 1;

        self.f = (self.f + 1) % DEG;
        self.rear = (self.rear + 1) % DEG;

        result
    }

    /// The C program's `flip()`.
    #[inline]
    pub fn flip(&mut self) -> bool {
        self.next_u32() & 1 == 1
    }
}

impl Default for GlibcRandom {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The first ten values of glibc `random()` after `srandom(1)` -- also the
    /// unseeded sequence, which is what the C program consumes.
    const GLIBC_SEED1: [u32; 10] = [
        1804289383, 846930886, 1681692777, 1714636915, 1957747793, 424238335, 719885386,
        1649760492, 596516649, 1189641421,
    ];

    #[test]
    fn matches_glibc_seed_1() {
        let mut rng = GlibcRandom::new();
        let got: Vec<u32> = (0..10).map(|_| rng.next_u32()).collect();
        assert_eq!(got, GLIBC_SEED1, "does not match glibc random() for seed 1");
    }

    #[test]
    fn explicit_seed_1_matches_unseeded() {
        let mut a = GlibcRandom::new();
        let mut b = GlibcRandom::with_seed(1);
        for _ in 0..100 {
            assert_eq!(a.next_u32(), b.next_u32());
        }
    }

    #[test]
    fn flip_reads_low_bit() {
        let mut rng = GlibcRandom::new();
        let expected: Vec<bool> = GLIBC_SEED1.iter().map(|v| v & 1 == 1).collect();
        let got: Vec<bool> = (0..10).map(|_| rng.flip()).collect();
        assert_eq!(got, expected);
    }

    /// Long-run check on the ring-pointer wraparound.
    ///
    /// These values are NOT an independent validation of the algorithm -- they
    /// come from a separate transcription of glibc's `random_r.c` (in Python),
    /// which shares this file's reading of glibc. They catch transcription slips
    /// in the Rust, nothing more. The only external anchor is the published
    /// seed-1 vector in `matches_glibc_seed_1`; final proof that this generator
    /// is the one the golden fixtures used arrives when the segmenter
    /// reproduces them (PLAN.md M3/M4).
    #[test]
    fn ring_wraparound_long_run() {
        let mut rng = GlibcRandom::new();
        let mut acc: u64 = 0;
        let mut at_9999 = 0;
        for i in 0..100_000u32 {
            let v = rng.next_u32();
            if i == 9999 {
                at_9999 = v;
            }
            acc = acc.wrapping_mul(1000003).wrapping_add(v as u64);
        }
        assert_eq!(at_9999, 1908609430, "draw 10000 diverged");
        assert_eq!(acc, 2636739091045304704, "100k-draw checksum diverged");
    }
}
