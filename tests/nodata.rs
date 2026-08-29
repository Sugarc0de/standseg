//! M5 gate: masking and nodata (water, non-treed area).
//!
//! Nodata does not get its own algorithm -- it is funnelled into the mask the
//! original already had. These tests pin the three properties that matter for
//! tree-stand segmentation: nodata pixels are excluded from the output, they
//! never leak into a centroid, and a stand can never grow across them.

use fast_segment::config::SegConfig;
use fast_segment::driver::{run, Observer, Phase};
use fast_segment::image::Image;
use fast_segment::region::RegionId;
use fast_segment::segment::Segmenter;

#[derive(Default)]
struct Capture {
    armap: Option<Vec<RegionId>>,
    centroids: Vec<f32>,
    nreg: usize,
}

impl Observer for Capture {
    fn on_map(&mut self, phase: Phase, _pass: usize, seg: &Segmenter) -> Result<(), String> {
        if phase == Phase::Auxiliary {
            self.armap = Some(seg.bands.rband.clone());
            self.centroids = seg.rl.ctr.clone();
            self.nreg = seg.nreg;
        }
        Ok(())
    }
}

fn img_from(nlines: usize, nsamps: usize, nbands: usize, data: Vec<u8>) -> Image {
    assert_eq!(data.len(), nlines * nsamps * nbands);
    let mut i = Image::new(nlines, nsamps, nbands);
    i.data = data;
    i
}

fn cfg() -> SegConfig {
    SegConfig {
        tols: vec![10.0],
        cm: 1.0,
        ..Default::default()
    }
}

fn segment(img: Image, mask: Option<&[u8]>) -> Capture {
    let mut c = Capture::default();
    run(img, &cfg(), mask, &mut c).expect("segmentation failed");
    c
}

/// A uniform field split by a masked stripe. Spectrally the two halves are
/// identical, so without the mask they merge into one region; with it they
/// cannot, because a stand may not grow across water.
fn striped() -> (Image, Vec<u8>) {
    let (nl, ns) = (9usize, 9usize);
    let data = vec![100u8; nl * ns];
    let mut mask = vec![1u8; nl * ns];
    for l in 0..nl {
        for s in 3..6 {
            mask[l * ns + s] = 0;
        }
    }
    (img_from(nl, ns, 1, data), mask)
}

#[test]
fn masked_pixels_are_region_zero() {
    let (img, mask) = striped();
    let out = segment(img, Some(&mask)).armap.expect("armap");
    for (p, &m) in mask.iter().enumerate() {
        if m == 0 {
            assert_eq!(out[p], 0, "masked pixel {p} was assigned a region");
        } else {
            assert_ne!(out[p], 0, "valid pixel {p} was left unassigned");
        }
    }
}

#[test]
fn a_region_never_grows_across_nodata() {
    let (img, mask) = striped();
    let out = segment(img, Some(&mask)).armap.expect("armap");

    let left = out[0 * 9 + 0];
    let right = out[0 * 9 + 8];
    assert_ne!(
        left, right,
        "a region spanned the nodata stripe -- pix_check_bounds_and_mask is not \
         sealing masked directions"
    );

    // And no single region may appear on both sides.
    for l in 0..9usize {
        for s in 0..3usize {
            for s2 in 6..9usize {
                assert_ne!(out[l * 9 + s], out[l * 9 + s2]);
            }
        }
    }
}

/// Without the mask the same image is one region, which is what makes the test
/// above discriminating rather than vacuous.
#[test]
fn unmasked_control_merges_into_one_region() {
    let (img, _) = striped();
    let c = segment(img, None);
    assert_eq!(c.nreg, 1, "control should merge to a single region");
}

/// Nodata values must never reach a centroid. Here nodata is 255 against a
/// uniform field of 10: any leakage pulls a boundary centroid off 10.0.
#[test]
fn nodata_never_contributes_to_a_centroid() {
    let (nl, ns) = (8usize, 8usize);
    let mut data = vec![10u8; nl * ns];
    let mut mask = vec![1u8; nl * ns];
    // A block of nodata in the corner, at a wildly different value.
    for l in 0..3 {
        for s in 0..3 {
            data[l * ns + s] = 255;
            mask[l * ns + s] = 0;
        }
    }
    let c = segment(img_from(nl, ns, 1, data), Some(&mask));
    for r in 1..=c.nreg {
        assert_eq!(
            c.centroids[r], 10.0,
            "region {r} centroid is {} -- nodata leaked in",
            c.centroids[r]
        );
    }
}

/// The derived-nodata path must agree with an equivalent explicit mask.
#[test]
fn derived_nodata_matches_explicit_mask() {
    let (nl, ns) = (8usize, 8usize);
    let mut data = vec![10u8; nl * ns];
    let mut mask = vec![1u8; nl * ns];
    for l in 2..5 {
        for s in 2..5 {
            data[l * ns + s] = 0;
            mask[l * ns + s] = 0;
        }
    }
    // Same rule the CLI applies for --nodata 0 with all-bands matching.
    let derived: Vec<u8> = data.iter().map(|&v| if v == 0 { 0 } else { 1 }).collect();
    assert_eq!(derived, mask);

    let a = segment(img_from(nl, ns, 1, data.clone()), Some(&mask))
        .armap
        .unwrap();
    let b = segment(img_from(nl, ns, 1, data), Some(&derived))
        .armap
        .unwrap();
    assert_eq!(a, b);
}

/// An entirely-nodata scene must terminate rather than spin or panic.
#[test]
fn all_nodata_terminates() {
    let (nl, ns) = (6usize, 6usize);
    let c = segment(
        img_from(nl, ns, 1, vec![0u8; nl * ns]),
        Some(&vec![0u8; nl * ns]),
    );
    assert_eq!(c.nreg, 0);
    assert!(c.armap.unwrap().iter().all(|&r| r == 0));
}

/// Multi-band nodata: the default rule is "all bands match".
#[test]
fn multiband_nodata_all_bands_rule() {
    let (nl, ns, nb) = (4usize, 4usize, 3usize);
    let mut data = vec![50u8; nl * ns * nb];
    // Pixel 0: all three bands are 0 -> nodata.
    for b in 0..nb {
        data[b] = 0;
    }
    // Pixel 1: only one band is 0 -> ordinary dark ground, still valid.
    data[nb] = 0;

    let all: Vec<u8> = (0..nl * ns)
        .map(|p| {
            if (0..nb).all(|b| data[p * nb + b] == 0) {
                0
            } else {
                1
            }
        })
        .collect();
    assert_eq!(all[0], 0, "pixel 0 should be nodata");
    assert_eq!(all[1], 1, "pixel 1 should NOT be nodata under the all-bands rule");

    let out = segment(img_from(nl, ns, nb, data), Some(&all)).armap.unwrap();
    assert_eq!(out[0], 0);
    assert_ne!(out[1], 0);
}
