//! Region list and the merge operation.
//!
//! Layout note: the C keeps an array of 12-byte `region` structs. We keep a
//! structure of arrays, because `reg_nnbr` and `merge_regions` touch bounding
//! boxes and centroids constantly but flags rarely. That is a locality change
//! only -- it cannot alter a result.
//!
//! Field widths follow PLAN.md section 4, trick 3: `u16` bounding-box
//! coordinates (the C's 32767 ceiling came from those being *signed*; nothing
//! ever stores a negative coordinate, so unsigned buys 65535 per axis for free).
//!
//! `npix` is `u32`, not the C's `unsigned short`. 65535 pixels is a 256 m square
//! at 1 m resolution -- smaller than plenty of real forest stands -- and the C
//! did not clamp there, it aborted the run. Widening costs 2 bytes per region
//! (226 MB at 15000^2, against a 5 GB peak) and removes a failure mode. See
//! PLAN.md section 12.2.

use crate::contig::Connectivity;
use crate::image::Sample;

pub const RF_ACTIVE: u8 = 1 << 0;
pub const RF_MERGE: u8 = 1 << 1;
pub const RF_SPECIAL: u8 = 1 << 2;

/// The real ceiling on a region's pixel count now: `u16` bounding-box
/// coordinates cap an image at 65536 x 65536, whose pixel count is exactly
/// `u32::MAX + 1`, so a region can hold at most `u32::MAX` pixels.
pub const MAX_REGION_PIXELS: u32 = u32::MAX;

pub type RegionId = u32;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct BBox {
    pub ulx: u16,
    pub uly: u16,
    pub lrx: u16,
    pub lry: u16,
}

/// One entry per region id. Index 0 is the artificial region holding masked and
/// nodata pixels; index `nreg + 1` is the scratch region `make_region_list` uses
/// to fold pixel pairs together.
pub struct RegionList {
    pub bbox: Vec<BBox>,
    pub npix: Vec<u32>,
    pub flags: Vec<u8>,
    /// Centroids, `nbands`-strided. f32 throughout, and repeatedly rounded --
    /// see PLAN.md section 3.2: this cannot be reconstructed from integer sums,
    /// which is why it dominates the memory budget.
    pub ctr: Vec<f32>,
    pub nbands: usize,
}

impl RegionList {
    /// `slots` includes index 0 and any scratch slot the caller wants.
    pub fn new(slots: usize, nbands: usize) -> Self {
        Self {
            bbox: vec![BBox::default(); slots],
            npix: vec![0; slots],
            flags: vec![0; slots],
            ctr: vec![0.0; slots * nbands],
            nbands,
        }
    }

    #[inline]
    pub fn is_active(&self, r: RegionId) -> bool {
        self.flags[r as usize] & RF_ACTIVE != 0
    }

    #[inline]
    pub fn ctr(&self, r: RegionId) -> &[f32] {
        let o = r as usize * self.nbands;
        &self.ctr[o..o + self.nbands]
    }

    /// Squared centroid distance, accumulated band by band in f32.
    ///
    /// The order and the width both matter: reassociating this or letting it
    /// contract into an FMA changes the low bits, which changes tie detection,
    /// which changes how many `flip()` draws are consumed. PLAN.md section 3.2.
    #[inline]
    pub fn dist2(&self, a: RegionId, b: RegionId) -> f32 {
        let (ca, cb) = (self.ctr(a), self.ctr(b));
        // zip rather than index: it drops a bounds check per band without
        // touching the accumulation order, which has to stay exactly as the C
        // has it (section 3.2).
        let mut dist2 = 0.0f32;
        for (x, y) in ca.iter().zip(cb.iter()) {
            let diff = *x - *y;
            dist2 += diff * diff;
        }
        dist2
    }

    /// The C's `region_from_pixel`: a fresh one-pixel region at (x, y).
    ///
    /// Generic over the sample width. For `u8` this is `pix[b] as f32`, exactly
    /// as the C has it; wider samples convert the same way and stay exact --
    /// every `u16` and `i16` is representable in f32.
    pub fn from_pixel<T: Sample>(&mut self, r: RegionId, x: u16, y: u16, pix: &[T]) {
        let i = r as usize;
        self.bbox[i] = BBox {
            ulx: x,
            uly: y,
            lrx: x,
            lry: y,
        };
        self.npix[i] = 1;
        self.flags[i] = RF_ACTIVE;
        let o = i * self.nbands;
        for b in 0..self.nbands {
            self.ctr[o + b] = pix[b].to_f32();
        }
    }
}

/// Merge `r2` into `r1`, deactivating `r2`.
///
/// Faithful to `region.c: merge_regions`, including the order of the three
/// stages: centroid, then contiguity over the *union* box while the region band
/// still carries both labels, then relabel over `r2`'s old box.
pub fn merge_regions(
    rl: &mut RegionList,
    rband: &mut [RegionId],
    cband: &mut [u8],
    nsamps: usize,
    conn: &Connectivity,
    offs: &[isize; 8],
    r1: RegionId,
    r2: RegionId,
) -> Result<(), String> {
    debug_assert_ne!(r1, r2);
    let (i1, i2) = (r1 as usize, r2 as usize);

    // Weighted centroid. The products go through int -> float promotion and the
    // divisor is computed in *int* before conversion; both mirror C's usual
    // arithmetic conversions.
    let n1 = rl.npix[i1];
    let n2 = rl.npix[i2];
    {
        let nb = rl.nbands;
        let (o1, o2) = (i1 * nb, i2 * nb);
        let denom = (n1 as u64 + n2 as u64) as f32;
        for b in 0..nb {
            let c1 = rl.ctr[o1 + b];
            let c2 = rl.ctr[o2 + b];
            rl.ctr[o1 + b] = (n1 as f32 * c1 + n2 as f32 * c2) / denom;
        }
    }

    let b2 = rl.bbox[i2];
    let nb1 = {
        let b1 = rl.bbox[i1];
        BBox {
            ulx: b1.ulx.min(b2.ulx),
            uly: b1.uly.min(b2.uly),
            lrx: b1.lrx.max(b2.lrx),
            lry: b1.lry.max(b2.lry),
        }
    };
    rl.bbox[i1] = nb1;

    let Some(mpix) = n1.checked_add(n2) else {
        return Err(format!(
            "merged region too large ({} pixels) from regions {r1} and {r2}; the \
             ceiling is {MAX_REGION_PIXELS}",
            n1 as u64 + n2 as u64
        ));
    };
    rl.npix[i1] = mpix;

    // Contiguity: any boundary pixel of either region that now abuts the other
    // becomes internal in that direction.
    //
    // Bits are only ever set here, and each direction tests its own bit, so
    // accumulating into a local and storing once is equivalent to the C's
    // repeated read-modify-write of `*Curmap`.
    let (ncdir, internal) = (conn.ncdir, conn.internal);
    let mut flags = [0u8; 8];
    flags[..ncdir].copy_from_slice(&conn.flags[..ncdir]);
    for l in nb1.uly as usize..=nb1.lry as usize {
        let row = l * nsamps;
        for p in row + nb1.ulx as usize..=row + nb1.lrx as usize {
            let rid = rband[p];
            if (rid != r1 && rid != r2) || cband[p] == internal {
                continue;
            }
            let mut map = cband[p];
            for d in 0..ncdir {
                if map & flags[d] == 0 {
                    // Safe without a bounds test: boundary pixels have their
                    // out-of-bounds directions marked contiguous by
                    // pix_check_bounds_and_mask, so a clear bit implies in-bounds.
                    let np = (p as isize + offs[d]) as usize;
                    let r = rband[np];
                    if r == r1 || r == r2 {
                        map |= flags[d];
                    }
                }
            }
            cband[p] = map;
        }
    }

    // Relabel r2's pixels.
    for l in b2.uly as usize..=b2.lry as usize {
        let row = l * nsamps;
        for v in &mut rband[row + b2.ulx as usize..=row + b2.lrx as usize] {
            if *v == r2 {
                *v = r1;
            }
        }
    }

    rl.flags[i2] = 0;
    Ok(())
}
