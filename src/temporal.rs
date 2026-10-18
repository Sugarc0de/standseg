//! Segmentation over a time series of images (Ye et al., paper 2).
//!
//! The two-stage segmenter of [`crate::driver`] answers "what are the segments
//! in this image". Over a series it is the wrong question to ask afresh each
//! time: re-segmenting every step independently moves boundaries that nothing on
//! the ground moved, so segments are not comparable between steps. This module
//! carries each step's segmentation forward and re-segments only where something
//! changed:
//!
//! | step | what happens |
//! |---|---|
//! | 0 | stage 1 on image 0, then stage 2 on layer 0 -- exactly the two-stage run |
//! | t >= 1 | stage 1 on image t **inside the change mask only**; paste those segments over the previous result ([`combine`]); stage 2 on layer t |
//!
//! The change mask is either given or derived from a categorical layer per step
//! (land cover, a burn-severity class, any disturbance product): a pixel is
//! re-segmented at step `t` when its class differs from step `t - 1`
//! ([`change_mask`]).
//!
//! Nothing here knows what the images depict. The oracle that defines it is the
//! Python of `forest-stand-update` (`run_our_algorithm.py` and its helpers),
//! where the series is annual Landsat proxies, the categorical layer is land
//! cover, and the second-stage layer is stem volume -- see `tests/TEMPORAL.md`.
//!
//! **Stage 1 sees nodata as data.** The C the oracle shelled out to never read
//! a nodata value, while [`crate::io::read_with_nodata`] honours one. The
//! series driver therefore takes its stage-1 mask from the change mask alone;
//! a caller that wants nodata masked as well has to say so.

use crate::config::SegConfig;
use crate::driver::{run_normal, Observer};
use crate::image::{Image, Samples};
use crate::stage2::{self, Stage2Config, Stage2Result};

/// Rescale every band to 1..=255, with nodata pixels set to 0.
///
/// This is the oracle's `normalize_image`, as `write_normalized_raster` applies
/// it to every input before segmenting. Per band: `low`/`high` are the band's
/// minimum and maximum over pixels that are not nodata, and each such pixel
/// becomes `trunc((x - low) / (high - low) * 254 + 1)`, computed in f64. The
/// truncation is the numpy cast back into the integer array; the `/`, `*` and
/// `+` are each one correctly rounded f64 operation in the same order, so the
/// result is bit-for-bit numpy's.
///
/// Nodata is per *band* (`img == nodata`, elementwise), not per pixel -- unlike
/// stage 1's all-bands rule.
///
/// Float input is refused: numpy would round the scaled value to f32 before the
/// cast to uint8, which this does not reproduce. A band that is constant over
/// its valid pixels has no range to rescale and is an error, as it was an
/// assertion failure in the oracle. With no nodata value at all the oracle
/// raised `NameError`; here every pixel is simply valid.
pub fn normalize(img: &Image, nodata: Option<f64>) -> Result<Image, String> {
    fn go<T: Copy + Into<f64>>(
        v: &[T],
        nbands: usize,
        nodata: Option<f64>,
    ) -> Result<Vec<u8>, String> {
        let mut out = vec![0u8; v.len()];
        for b in 0..nbands {
            let valid = |x: f64| nodata != Some(x);
            let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
            for s in v.iter().skip(b).step_by(nbands) {
                let x: f64 = (*s).into();
                if valid(x) {
                    lo = lo.min(x);
                    hi = hi.max(x);
                }
            }
            if !(hi > lo) {
                return Err(format!(
                    "band {} has no range to rescale (every valid pixel is {lo})",
                    b + 1
                ));
            }
            let span = hi - lo;
            for (i, s) in v.iter().enumerate().skip(b).step_by(nbands) {
                let x: f64 = (*s).into();
                if valid(x) {
                    // The integer subtraction is exact in the source type, and
                    // exact again in f64: every operand is below 2^16.
                    out[i] = ((x - lo) / span * 254.0 + 1.0) as u8;
                }
            }
        }
        Ok(out)
    }
    let data = match &img.data {
        Samples::U8(v) => go(v, img.nbands, nodata)?,
        Samples::U16(v) => go(v, img.nbands, nodata)?,
        Samples::I16(v) => go(v, img.nbands, nodata)?,
        Samples::F32(_) => {
            return Err("rescaling to 1..255 is defined for integer images only; this one is \
                 32-bit float"
                .into())
        }
    };
    let mut out = Image::from_samples(img.nlines, img.nsamps, img.nbands, Samples::U8(data));
    out.geo = img.geo.clone();
    Ok(out)
}

/// 1 where the categorical layer differs between two steps, 0 where it does not.
///
/// The oracle compares band 1 only (`read(1)`); a multi-band class layer here
/// counts a pixel as changed when any band differs, which is the same thing for
/// the single-band layers it was written for. Values are compared raw, so a
/// pixel that is nodata at both steps is unchanged.
pub fn change_mask(prev: &Image, cur: &Image) -> Result<Vec<u8>, String> {
    if (prev.nlines, prev.nsamps, prev.nbands) != (cur.nlines, cur.nsamps, cur.nbands) {
        return Err(format!(
            "class layers differ in shape: {}x{}x{} then {}x{}x{}",
            prev.nlines, prev.nsamps, prev.nbands, cur.nlines, cur.nsamps, cur.nbands
        ));
    }
    let nb = cur.nbands;
    fn diff<T: PartialEq>(a: &[T], b: &[T], nb: usize) -> Vec<u8> {
        a.chunks_exact(nb)
            .zip(b.chunks_exact(nb))
            .map(|(x, y)| u8::from(x != y))
            .collect()
    }
    Ok(match (&prev.data, &cur.data) {
        (Samples::U8(a), Samples::U8(b)) => diff(a, b, nb),
        (Samples::U16(a), Samples::U16(b)) => diff(a, b, nb),
        (Samples::I16(a), Samples::I16(b)) => diff(a, b, nb),
        // Bitwise, so NaN == NaN: an unchanged nodata NaN is unchanged.
        (Samples::F32(a), Samples::F32(b)) => a
            .chunks_exact(nb)
            .zip(b.chunks_exact(nb))
            .map(|(x, y)| u8::from(x.iter().zip(y).any(|(p, q)| p.to_bits() != q.to_bits())))
            .collect(),
        _ => {
            return Err(format!(
                "class layers differ in sample type: {} then {}",
                prev.data.kind(),
                cur.data.kind()
            ))
        }
    })
}

/// Paste `update` over `prev` wherever `update` is nonzero, then give every
/// 4-connected piece of every id its own id.
///
/// The oracle's `combine_two_region_map_with_no_change_mask` (as fixed in
/// `981892d`). It shifts the update's ids past `prev`'s maximum, pastes, then
/// visits the distinct ids in ascending order and, within each, labels the
/// 4-connected components with `scipy.ndimage.label` -- which numbers them in
/// raster order of each component's first pixel -- and finally renumbers
/// everything 1..n in that order. So the output id of a piece is its rank under
/// (pasted id, raster position of its first pixel), plus one; 0 stays 0.
///
/// Why split pieces at all: an old segment the change cut in two would
/// otherwise be one id in two places, and stage 2 would treat the halves as one
/// region with one centroid.
///
/// One departure, on input the oracle cannot have met: it renumbers with
/// `enumerate(np.unique(map))`, which gives 0 to whatever the smallest value is,
/// so a combined map with *no* 0 pixel would turn a real segment into masked
/// area. Here 0 is only ever 0. Our oracle run asserts the situation never
/// arises on the data pinned in `tests/temporal/`.
///
/// Returns the map and the number of nonzero ids in it.
pub fn combine(
    prev: &[u32],
    update: &[u32],
    nlines: usize,
    nsamps: usize,
) -> Result<(Vec<u32>, usize), String> {
    let n = nlines * nsamps;
    if prev.len() != n || update.len() != n {
        return Err("combine: maps are not the same size".into());
    }
    let max_prev = prev.iter().copied().max().unwrap_or(0) as u64;
    let key = |p: usize| -> u64 {
        match update[p] {
            0 => prev[p] as u64,
            u => u as u64 + max_prev,
        }
    };

    // Components in raster order of their first pixel: scanning in raster
    // order, the first unlabelled pixel met *is* its component's first pixel.
    const NONE: u32 = u32::MAX;
    let mut comp = vec![NONE; n];
    let mut first_key: Vec<u64> = Vec::new();
    let mut stack: Vec<usize> = Vec::new();
    for start in 0..n {
        let k = key(start);
        if k == 0 || comp[start] != NONE {
            continue;
        }
        let c = first_key.len() as u32;
        if c == NONE {
            return Err("combine: more than 2^32 - 1 pieces".into());
        }
        first_key.push(k);
        comp[start] = c;
        stack.push(start);
        while let Some(p) = stack.pop() {
            let (l, s) = (p / nsamps, p % nsamps);
            let mut visit = |q: usize| {
                if comp[q] == NONE && key(q) == k {
                    comp[q] = c;
                    stack.push(q);
                }
            };
            if l > 0 {
                visit(p - nsamps);
            }
            if s + 1 < nsamps {
                visit(p + 1);
            }
            if l + 1 < nlines {
                visit(p + nsamps);
            }
            if s > 0 {
                visit(p - 1);
            }
        }
    }

    // Rank by (key, discovery order). A stable sort on key keeps discovery
    // order within a key.
    let mut order: Vec<u32> = (0..first_key.len() as u32).collect();
    order.sort_by_key(|&c| first_key[c as usize]);
    let mut new_id = vec![0u32; first_key.len()];
    for (rank, &c) in order.iter().enumerate() {
        new_id[c as usize] = rank as u32 + 1;
    }
    let out = comp
        .iter()
        .map(|&c| if c == NONE { 0 } else { new_id[c as usize] })
        .collect();
    Ok((out, first_key.len()))
}

/// Bytes per id the C's writer would choose for `nreg` regions (the bit width
/// of the count, rounded up to 1, 2 or 4 bytes).
pub fn nbytes_for(nreg: usize) -> usize {
    match usize::BITS - nreg.leading_zeros() {
        0..=8 => 1,
        9..=16 => 2,
        _ => 4,
    }
}

/// What one step produced, for the caller to write and report.
pub struct StepResult {
    /// The step-t stage-1 map: the whole image at step 0, only the changed area
    /// (0 elsewhere) after.
    pub stage1: Vec<u32>,
    pub stage1_passes: usize,
    pub stage1_nbytes: usize,
    /// What stage 2 started from: the stage-1 map at step 0, the [`combine`]d
    /// map after. `None` at step 0, where it would only repeat `stage1`.
    pub combined: Option<Vec<u32>>,
    /// The step's result.
    pub map: Vec<u32>,
    pub stage2: Stage2Result,
    /// Width every map in the series is written at; see [`Series`].
    pub nbytes: usize,
    /// Pixels re-segmented this step (all of them at step 0).
    pub changed: usize,
}

/// The state carried from one step to the next: the previous result and the
/// width it is written at.
///
/// **Width.** The oracle writes each step's maps with the previous map's
/// rasterio profile, so the id width is fixed by the first stage-1 map (the C's
/// rule, from its region count) and inherited by every step after -- even when
/// a later map would fit in fewer bytes. Kept, so payloads compare byte for
/// byte. Unlike the oracle, a map whose ids outgrow the inherited width is
/// widened rather than wrapped.
pub struct Series<'c> {
    pub cfg: &'c SegConfig,
    pub s2: Stage2Config,
    pub nlines: usize,
    pub nsamps: usize,
    prev: Option<Vec<u32>>,
    nbytes: usize,
}

impl<'c> Series<'c> {
    pub fn new(cfg: &'c SegConfig, s2: Stage2Config, nlines: usize, nsamps: usize) -> Self {
        Self {
            cfg,
            s2,
            nlines,
            nsamps,
            prev: None,
            nbytes: 0,
        }
    }

    /// Run the next step. `change` must be `None` for the first step and given
    /// for every step after it (1 = re-segment this pixel).
    pub fn step(
        &mut self,
        image: Image,
        change: Option<&[u8]>,
        layer: &Image,
        obs: &mut dyn Observer,
    ) -> Result<StepResult, String> {
        let (nl, ns) = (self.nlines, self.nsamps);
        for (what, i) in [("image", &image), ("second-stage layer", layer)] {
            if (i.nlines, i.nsamps) != (nl, ns) {
                return Err(format!(
                    "{what} is {}x{}, but the series is {nl}x{ns}",
                    i.nlines, i.nsamps
                ));
            }
        }
        let npix = nl * ns;
        match (&self.prev, change) {
            (None, Some(_)) => return Err("the first step has nothing to change from".into()),
            (Some(_), None) => return Err("every step after the first needs a change mask".into()),
            (_, Some(m)) if m.len() != npix => {
                return Err(format!("change mask is {} pixels, expected {npix}", m.len()))
            }
            _ => {}
        }
        let changed = change.map_or(npix, |m| m.iter().filter(|&&v| v != 0).count());

        let (seg, stage1_passes) = run_normal(image, self.cfg, change, obs)?;
        let stage1_nbytes = seg.region_map_nbytes();
        let stage1 = seg.bands.rband[..npix].to_vec();
        drop(seg);

        let (mut map, combined) = match self.prev.take() {
            None => {
                self.nbytes = stage1_nbytes;
                (stage1.clone(), None)
            }
            Some(prev) => {
                let (c, nreg) = combine(&prev, &stage1, nl, ns)?;
                self.nbytes = self.nbytes.max(nbytes_for(nreg));
                (c.clone(), Some(c))
            }
        };
        let stage2 = stage2::run(&mut map, nl, ns, layer, &self.s2)?;
        self.prev = Some(map.clone());
        Ok(StepResult {
            stage1,
            stage1_passes,
            stage1_nbytes,
            combined,
            map,
            stage2,
            nbytes: self.nbytes,
            changed,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn combine_shifts_pastes_and_splits() {
        // prev: id 5 is cut in two by the update; id 7 untouched; 0 stays 0.
        #[rustfmt::skip]
        let prev = [
            5, 5, 5, 7,
            5, 5, 5, 7,
            0, 0, 7, 7,
        ];
        #[rustfmt::skip]
        let update = [
            0, 1, 0, 0,
            0, 1, 0, 0,
            0, 2, 0, 0,
        ];
        let (out, n) = combine(&prev, &update, 3, 4).unwrap();
        // keys: 5 (two pieces: left col, right col), 7, then 1+7=8, 2+7=9.
        // ranks: 5/left=1, 5/right=2, 7=3, 8=4, 9=5.
        #[rustfmt::skip]
        let want = [
            1, 4, 2, 3,
            1, 4, 2, 3,
            0, 5, 3, 3,
        ];
        assert_eq!(out, want);
        assert_eq!(n, 5);
    }

    #[test]
    fn combine_pieces_follow_first_pixel_not_id() {
        // Id 9 appears first in raster order but sorts after id 3.
        let prev = [9, 0, 3, 0, 9];
        let update = [0; 5];
        let (out, _) = combine(&prev, &update, 1, 5).unwrap();
        assert_eq!(out, [2, 0, 1, 0, 3]);
    }

    #[test]
    fn combine_keeps_zero_when_there_is_no_zero() {
        // The oracle would hand 0 to the smallest id here; we do not.
        let (out, n) = combine(&[4, 4, 6], &[0, 0, 0], 1, 3).unwrap();
        assert_eq!(out, [1, 1, 2]);
        assert_eq!(n, 2);
    }

    #[test]
    fn normalize_matches_numpy_rule() {
        // band 0: valid 10..=20, nodata 65535; band 1: valid 0..=3.
        let v: Vec<u16> = vec![10, 0, 65535, 3, 15, 1, 20, 2];
        let img = Image::from_samples(1, 4, 2, Samples::U16(v));
        let out = normalize(&img, Some(65535.0)).unwrap();
        // (15-10)/10*254+1 = 128.0 ; (1-0)/3*254+1 = 85.666 -> 85 ; 2/3 -> 170.33 -> 170
        assert_eq!(out.data, Samples::U8(vec![1, 1, 0, 255, 128, 85, 255, 170]));
    }

    #[test]
    fn change_mask_is_elementwise() {
        let a = Image::from_samples(1, 3, 1, Samples::U8(vec![210, 255, 80]));
        let b = Image::from_samples(1, 3, 1, Samples::U8(vec![20, 255, 80]));
        assert_eq!(change_mask(&a, &b).unwrap(), [1, 0, 0]);
    }

    #[test]
    fn nbytes_rule() {
        assert_eq!(nbytes_for(0), 1);
        assert_eq!(nbytes_for(255), 1);
        assert_eq!(nbytes_for(256), 2);
        assert_eq!(nbytes_for(65535), 2);
        assert_eq!(nbytes_for(65536), 4);
    }
}
