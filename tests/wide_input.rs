//! Wide-sample input: 16-bit unsigned and signed.
//!
//! The 1992 program was uint8-only. Landsat 8/9 and Sentinel-2 are 12-bit data
//! in a 16-bit container, so rescaling to 8 bits before segmenting throws away
//! radiometry and changes the answer. These tests pin three things:
//!
//! 1. widening a byte image to 16 bits, values unchanged, gives *bit-identical*
//!    region maps -- so the wide path is a generalisation, not a second
//!    algorithm with its own behaviour;
//! 2. real 16-bit imagery segments, and does not agree with its own 8-bit
//!    rescaling (which is the reason to support it at all);
//! 3. a negative nodata sentinel works, because int16 is how Landsat
//!    Collection 2 surface reflectance actually ships and -9999 is its fill.

use std::path::{Path, PathBuf};

use fast_segment::config::SegConfig;
use fast_segment::driver::{run, Observer, Phase};
use fast_segment::image::{Image, Samples};
use fast_segment::region::RegionId;
use fast_segment::segment::Segmenter;

fn golden(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden").join(rel)
}

#[derive(Default)]
struct Capture {
    rmap: Option<Vec<RegionId>>,
    armap: Option<Vec<RegionId>>,
    apass: usize,
    nreg: usize,
}

impl Observer for Capture {
    fn on_map(&mut self, phase: Phase, pass: usize, seg: &Segmenter) -> Result<(), String> {
        match phase {
            Phase::Normal => self.rmap = Some(seg.bands.rband.clone()),
            Phase::Auxiliary => {
                self.armap = Some(seg.bands.rband.clone());
                self.apass = pass;
                self.nreg = seg.nreg;
            }
        }
        Ok(())
    }
}

fn segment(img: Image, cfg: &SegConfig, mask: Option<&[u8]>) -> Capture {
    let mut cap = Capture::default();
    run(img, cfg, mask, &mut cap).expect("segmentation failed");
    cap
}

fn case1_config() -> SegConfig {
    SegConfig {
        tols: vec![10.0],
        cm: 0.1,
        ..Default::default()
    }
    .with_n(&[15, 15, 100, 2500, 2500])
    .unwrap()
}

/// Widen every sample u8 -> u16 without changing its value.
fn widen(img: &Image) -> Image {
    let v: Vec<u16> = img
        .data
        .as_u8()
        .expect("source must be 8-bit")
        .iter()
        .map(|&b| u16::from(b))
        .collect();
    Image::from_samples(img.nlines, img.nsamps, img.nbands, Samples::U16(v))
}

/// The load-bearing test. Case 1 is the fixture whose byte-exact reproduction
/// is the project's definition of done; run it again with every sample stored
/// in 16 bits instead of 8 and the region maps must be identical, id for id.
/// If the wide path drifted -- a different accumulation width, a lost tie, one
/// extra `flip()` draw -- this is where it shows.
#[test]
fn widening_case1_to_u16_changes_nothing() {
    let narrow = fast_segment::io::read(&golden("misc/temp_byte_bip")).expect("read");
    let wide = widen(&narrow);

    let cfg = case1_config();
    let a = segment(narrow, &cfg, None);
    let b = segment(wide, &cfg, None);

    assert_eq!(
        a.rmap, b.rmap,
        "rmap differs between the 8-bit and 16-bit paths on identical values"
    );
    assert_eq!(
        a.armap, b.armap,
        "armap differs between the 8-bit and 16-bit paths on identical values"
    );
    assert_eq!(a.apass, b.apass, "auxiliary pass count differs");
    assert_eq!(a.nreg, b.nreg, "final region count differs");
}

/// Same again for int16, which additionally exercises a signed sample type.
#[test]
fn widening_case1_to_i16_changes_nothing() {
    let narrow = fast_segment::io::read(&golden("misc/temp_byte_bip")).expect("read");
    let v: Vec<i16> = narrow
        .data
        .as_u8()
        .unwrap()
        .iter()
        .map(|&b| i16::from(b))
        .collect();
    let wide = Image::from_samples(
        narrow.nlines,
        narrow.nsamps,
        narrow.nbands,
        Samples::I16(v),
    );

    let cfg = case1_config();
    let a = segment(narrow, &cfg, None);
    let b = segment(wide, &cfg, None);
    assert_eq!(a.armap, b.armap, "armap differs between the u8 and i16 paths");
}

/// The real 16-bit Landsat stack: it segments, every pixel gets a label, and
/// the answer is *not* the answer you get from the 8-bit rescaling of the same
/// scene. That difference is the reason for this whole change.
#[test]
fn int16_landsat_stack_segments_and_differs_from_its_8bit_rescaling() {
    let dir = golden("LC80220492014083LGN00/input");
    let wide = fast_segment::io::read(&dir.join("LC80220492014083LGN00_stack")).expect("read");
    let byte = fast_segment::io::read(&dir.join("LC80220492014083LGN00_stack.ipw")).expect("read");
    assert_eq!(wide.data.kind(), "16-bit signed");

    // Tolerance is in DN. The stack runs 0..8990 where the .ipw runs 0..255, so
    // a comparable tolerance is scaled by roughly the same factor.
    let wide_cfg = SegConfig {
        tols: vec![350.0],
        cm: 0.1,
        ..Default::default()
    }
    .with_n(&[15, 15, 100, 2500, 2500])
    .unwrap();

    let w = segment(wide, &wide_cfg, None);
    let b = segment(byte, &case1_config(), None);

    let wmap = w.armap.expect("no armap from the 16-bit run");
    let bmap = b.armap.expect("no armap from the 8-bit run");
    assert_eq!(wmap.len(), 250 * 250);
    assert!(w.nreg > 0, "16-bit run produced no regions");
    assert!(
        wmap.iter().all(|&r| r > 0),
        "every pixel should be labelled when nothing is masked"
    );
    assert_ne!(
        wmap, bmap,
        "segmenting 16-bit reflectance gave the same map as its 8-bit rescaling, \
         which would mean the extra radiometry was being discarded"
    );
}

/// -9999 is the Landsat Collection 2 fill value. A negative nodata sentinel has
/// to reach the mask, and masked pixels have to come back as region 0.
#[test]
fn negative_nodata_sentinel_masks_int16_fill() {
    let (nl, ns, nb) = (8usize, 8usize, 2usize);
    let mut v = vec![0i16; nl * ns * nb];
    for p in 0..nl * ns {
        let fill = p % ns == 0;
        for b in 0..nb {
            v[p * nb + b] = if fill { -9999 } else { 1000 + (p as i16 % 3) };
        }
    }
    let img = Image::from_samples(nl, ns, nb, Samples::I16(v));

    let mut mask = vec![1u8; nl * ns];
    img.apply_nodata(-9999, false, &mut mask);
    assert_eq!(
        mask.iter().filter(|&&m| m == 0).count(),
        nl,
        "one fill pixel per row should be masked"
    );

    let cfg = SegConfig {
        tols: vec![50.0],
        cm: 1.0,
        ..Default::default()
    };
    let cap = segment(img, &cfg, Some(&mask));
    let armap = cap.armap.expect("no armap");
    for p in 0..nl * ns {
        if mask[p] == 0 {
            assert_eq!(armap[p], 0, "masked fill pixel {p} is not region 0");
        } else {
            assert_ne!(armap[p], 0, "valid pixel {p} was left unlabelled");
        }
    }
}
