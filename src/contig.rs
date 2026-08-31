//! The contiguity band: one byte per pixel, one bit per direction.
//!
//! This is a real space optimisation, not a C workaround -- 8x smaller than a
//! byte per direction, and 225 MB rather than 1.8 GB on a 15000^2 image. What we
//! do *not* keep is the C's habit of silently changing what the byte means three
//! times over the run (`pixel.c`'s 30-line comment). PLAN.md section 4, trick 2.
//!
//! The three meanings, in order:
//!   1. after `pix_nnbr`  -- "this neighbour is at my minimum distance"
//!   2. after `pix_merge` -- "this is the one neighbour I merged with" (<=1 bit)
//!   3. after `pix_check_bounds_and_mask` and thereafter --
//!      "this neighbour is out of bounds, masked, or in my region", so a clear
//!      bit means "in bounds, and belongs to some *other* region".

pub const DF_N: u8 = 1 << 0;
pub const DF_NE: u8 = 1 << 1;
pub const DF_E: u8 = 1 << 2;
pub const DF_SE: u8 = 1 << 3;
pub const DF_S: u8 = 1 << 4;
pub const DF_SW: u8 = 1 << 5;
pub const DF_W: u8 = 1 << 6;
pub const DF_NW: u8 = 1 << 7;

pub const CINTERNAL4: u8 = DF_N | DF_E | DF_S | DF_W;
pub const CINTERNAL8: u8 = 0xff;
pub const CMONO: u8 = 0;
pub const CCLEAR: u8 = CMONO;

/// Boundary masks. Applied whole regardless of connectivity, exactly as the C
/// does -- in 4-way mode the diagonal bits are inert for `has_contig` but still
/// change the `== Cinternal` test, and that difference is load-bearing.
pub const N_EDGE: u8 = DF_NW | DF_N | DF_NE;
pub const E_EDGE: u8 = DF_NE | DF_E | DF_SE;
pub const S_EDGE: u8 = DF_SE | DF_S | DF_SW;
pub const W_EDGE: u8 = DF_SW | DF_W | DF_NW;

/// 4- or 8-way connectivity: direction offsets and their bit flags.
#[derive(Debug, Clone, Copy)]
pub struct Connectivity {
    pub ncdir: usize,
    pub deltas: &'static [(i32, i32)],
    pub flags: &'static [u8],
    pub internal: u8,
}

/// `(dx, dy)`, matching the C's `pcoord { short x; short y; }` order.
static CD4_DELTA: [(i32, i32); 4] = [(0, -1), (1, 0), (0, 1), (-1, 0)];
static CD4_FLAG: [u8; 4] = [DF_N, DF_E, DF_S, DF_W];

static CD8_DELTA: [(i32, i32); 8] = [
    (0, -1),
    (1, -1),
    (1, 0),
    (1, 1),
    (0, 1),
    (-1, 1),
    (-1, 0),
    (-1, -1),
];
static CD8_FLAG: [u8; 8] = [DF_N, DF_NE, DF_E, DF_SE, DF_S, DF_SW, DF_W, DF_NW];

pub const FOUR_WAY: Connectivity = Connectivity {
    ncdir: 4,
    deltas: &CD4_DELTA,
    flags: &CD4_FLAG,
    internal: CINTERNAL4,
};

pub const EIGHT_WAY: Connectivity = Connectivity {
    ncdir: 8,
    deltas: &CD8_DELTA,
    flags: &CD8_FLAG,
    internal: CINTERNAL8,
};

impl Connectivity {
    #[inline]
    pub fn has(&self, map: u8, d: usize) -> bool {
        map & self.flags[d] != 0
    }

    #[inline]
    pub fn set(&self, map: &mut u8, d: usize) {
        *map |= self.flags[d];
    }

    /// Flat neighbour offsets for a given row stride, so the hot loops can do
    /// one add instead of recomputing `dy * nsamps + dx` per pixel per direction.
    pub fn offsets(&self, nsamps: usize) -> [isize; 8] {
        let mut o = [0isize; 8];
        for (d, &(dx, dy)) in self.deltas.iter().enumerate().take(self.ncdir) {
            o[d] = dy as isize * nsamps as isize + dx as isize;
        }
        o
    }

    /// The C's `dir_reverse`: `(d + Ncdir/2) % Ncdir`.
    #[inline]
    pub fn reverse(&self, d: usize) -> usize {
        (d + self.ncdir / 2) % self.ncdir
    }
}
