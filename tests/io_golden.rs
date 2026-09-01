//! M1 gate: the readers agree with each other and with the fixtures on disk.

use std::path::{Path, PathBuf};

fn golden(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(rel)
}

/// The Case 1 input exists in two containers. If the ENVI and IPW readers do
/// not produce identical pixel buffers, one of them is wrong -- and the
/// segmentation would silently diverge depending on which file you passed.
#[test]
fn envi_and_ipw_readers_agree_on_case1() {
    let envi = standseg::io::read(&golden("misc/temp_byte_bip")).expect("read ENVI");
    let ipw = standseg::io::read(&golden("test_3456/input/test_3456.bip.ipw")).expect("read IPW");

    assert_eq!((envi.nlines, envi.nsamps, envi.nbands), (250, 250, 4));
    assert_eq!((ipw.nlines, ipw.nsamps, ipw.nbands), (250, 250, 4));
    assert_eq!(
        envi.data, ipw.data,
        "ENVI and IPW readers disagree on Case 1 pixels"
    );
}

#[test]
fn reads_case2_ipw() {
    let img = standseg::io::read(&golden(
        "LC80220492014083LGN00/input/LC80220492014083LGN00_stack.ipw",
    ))
    .expect("read Case 2 IPW");
    assert_eq!((img.nlines, img.nsamps, img.nbands), (250, 250, 8));
    assert_eq!(img.data.len(), 250 * 250 * 8);
}

/// Case 2's ENVI `_stack` is int16 -- the real Landsat reflectance, DN 0..8990.
/// The original program rejected it outright, which is why the golden was made
/// from the 8-bit `.ipw` rescaling instead. We can now read it, so the thing
/// worth pinning is no longer "it errors" but "it is not the same picture":
/// substituting it for the `.ipw` would silently produce a different answer.
#[test]
fn int16_stack_is_not_the_case2_input() {
    let dir = golden("LC80220492014083LGN00/input");
    let wide = standseg::io::read(&dir.join("LC80220492014083LGN00_stack"))
        .expect("int16 ENVI should now be readable");
    let ipw =
        standseg::io::read(&dir.join("LC80220492014083LGN00_stack.ipw")).expect("read Case 2 IPW");

    assert_eq!((wide.nlines, wide.nsamps, wide.nbands), (250, 250, 8));
    let w = wide.data.as_i16().expect("ENVI data type 2 reads as i16");
    let b = ipw.data.as_u8().expect("IPW band is 1 byte");
    assert_eq!(w.len(), b.len());
    assert!(
        w.iter().any(|&x| x > 255),
        "the _stack should hold 16-bit reflectance, not rescaled bytes"
    );
    let agree = w
        .iter()
        .zip(b)
        .filter(|(a, b)| i64::from(**a) == i64::from(**b))
        .count();
    assert!(
        agree * 2 < w.len(),
        "the int16 _stack unexpectedly resembles the .ipw -- re-read tests/GOLDEN.md \
         before using either as the Case 2 input"
    );
}

/// `test_3456.bip` shares a name with the real input but holds different data.
/// Pinning this stops anyone from quietly "fixing" a test by swapping inputs.
#[test]
fn test_3456_bip_is_not_the_case1_input() {
    let decoy = standseg::io::read(&golden("test_3456/input/test_3456.bip")).expect("read");
    let real = standseg::io::read(&golden("misc/temp_byte_bip")).expect("read");
    assert_ne!(
        decoy.data, real.data,
        "test_3456.bip unexpectedly matches the real input -- re-read tests/GOLDEN.md"
    );
}

/// Round-trip the region-map writer against a real golden output.
#[test]
fn region_map_writer_reproduces_golden_bytes() {
    // Load the golden armap payload as u16 region ids, write them back out, and
    // require the bytes to match. This validates the writer against the exact
    // artifact we have to reproduce at M4.
    let payload = std::fs::read(golden("test_3456/expected/proof/regmap.armap.58")).unwrap();
    assert_eq!(payload.len(), 125_000);
    let rband: Vec<u32> = payload
        .as_chunks::<2>()
        .0
        .iter()
        .map(|c| u16::from_le_bytes(*c) as u32)
        .collect();

    let dir = std::env::temp_dir().join("standseg_io_test");
    std::fs::create_dir_all(&dir).unwrap();
    let out = dir.join("roundtrip.armap.58");
    standseg::io::envi::write_region_map(
        &out,
        &rband,
        250,
        250,
        2,
        &standseg::image::GeoRef::default(),
        false,
        &standseg::io::Provenance::default(),
    )
    .expect("write");

    let got = std::fs::read(&out).unwrap();
    assert_eq!(
        got, payload,
        "region map writer does not reproduce golden bytes"
    );
    let _ = std::fs::remove_dir_all(&dir);
}
