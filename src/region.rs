//! Region list and the merge operation.
//!
//! Layout note: the C keeps an array of 12-byte `region` structs. We keep a
//! structure of arrays, because `reg_nnbr` and `merge_regions` touch bounding
//! boxes and centroids constantly but flags rarely. That is a locality change
//! only -- it cannot alter a result.
//!
//! Field widths follow PLAN.md section 4, trick 3: `u16` bounding-box
//! coordinates (the C's 32767 ceiling came from those being *signed*; nothing
//! ever stores a negative coordinate, so unsigned buys 65535 per axis for free)
//! and `u16` npix (already implied by the CLI validating `nabsmax <= 65535`).

use crate::contig::Connectivity;

pub const RF_ACTIVE: u8 = 1 << 0;
pub const RF_MERGE: u8 = 1 << 1;
pub const RF_SPECIAL: u8 = 1 << 2;

pub const MAX_USHORT: u32 = 65535;

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
    pub npix: Vec<u16>,
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
        let mut dist2 = 0.0f32;
        for i in 0..self.nbands {
            let diff = ca[i] - cb[i];
            dist2 += diff * diff;
        }
        dist2
    }

    /// The C's `region_from_pixel`: a fresh one-pixel region at (x, y).
    pub fn from_pixel(&mut self, r: RegionId, x: u16, y: u16, pix: &[u8]) {
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
            self.ctr[o + b] = pix[b] as f32;
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
    r1: RegionId,
    r2: RegionId,
) -> Result<(), String> {
    debug_assert_ne!(r1, r2);
    let (i1, i2) = (r1 as usize, r2 as usize);

    // Weighted centroid. The products go through int -> float promotion and the
    // divisor is computed in *int* before conversion; both mirror C's usual
    // arithmetic conversions.
    let n1 = rl.npix[i1] as u32;
    let n2 = rl.npix[i2] as u32;
    {
        let nb = rl.nbands;
        let (o1, o2) = (i1 * nb, i2 * nb);
        let denom = (n1 + n2) as f32;
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

    let mpix = n1 + n2;
    if mpix > MAX_USHORT {
        return Err(format!(
            "merged region too large ({mpix} pixels) from regions {r1} and {r2}"
        ));
    }
    rl.npix[i1] = mpix as u16;

    // Contiguity: any boundary pixel of either region that now abuts the other
    // becomes internal in that direction.
    for l in nb1.uly as usize..=nb1.lry as usize {
        for s in nb1.ulx as usize..=nb1.lrx as usize {
            let p = l * nsamps + s;
            let rid = rband[p];
            if rid != r1 && rid != r2 {
                continue;
            }
            if cband[p] == conn.internal {
                continue;
            }
            for d in 0..conn.ncdir {
                if !conn.has(cband[p], d) {
                    let (dx, dy) = conn.deltas[d];
                    // Safe without a bounds test: boundary pixels have their
                    // out-of-bounds directions marked contiguous by
                    // pix_check_bounds_and_mask, so a clear bit implies in-bounds.
                    let np = (l as i32 + dy) as usize * nsamps + (s as i32 + dx) as usize;
                    let r = rband[np];
                    if r == r1 || r == r2 {
                        conn.set(&mut cband[p], d);
                    }
                }
            }
        }
    }

    // Relabel r2's pixels.
    for l in b2.uly as usize..=b2.lry as usize {
        for s in b2.ulx as usize..=b2.lrx as usize {
            let p = l * nsamps + s;
            if rband[p] == r2 {
                rband[p] = r1;
            }
        }
    }

    rl.flags[i2] = 0;
    Ok(())
}
