//! M1b gate: TIFF and PNG readers.
//!
//! Rather than assert on pixel buffers, these re-encode the real Case 1 input
//! into each container and require the full segmentation to still land on the
//! golden bytes. A reader that transposed, padded or reordered anything would
//! fail loudly.

use std::path::{Path, PathBuf};

use fast_segment::config::SegConfig;
use fast_segment::driver::{run, Observer, Phase};
use fast_segment::segment::Segmenter;

fn golden(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden").join(rel)
}

fn tmp(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join("fast_segment_tiffpng");
    std::fs::create_dir_all(&d).unwrap();
    d.join(name)
}

fn case1() -> fast_segment::image::Image {
    fast_segment::io::read(&golden("misc/temp_byte_bip")).expect("read Case 1")
}

#[derive(Default)]
struct Capture {
    rmap: Option<(usize, Vec<u8>)>,
    armap: Option<(usize, Vec<u8>)>,
}

impl Observer for Capture {
    fn on_map(&mut self, phase: Phase, pass: usize, seg: &Segmenter) -> Result<(), String> {
        let nb = seg.region_map_nbytes();
        let mut out = Vec::new();
        for &v in &seg.bands.rband {
            out.extend_from_slice(&v.to_le_bytes()[..nb]);
        }
        match phase {
            Phase::Normal => self.rmap = Some((pass, out)),
            Phase::Auxiliary => self.armap = Some((pass, out)),
        }
        Ok(())
    }
}

/// Segment whatever is at `path` and require the Case 1 golden result.
fn assert_reproduces_golden(path: &Path, what: &str) {
    let img = fast_segment::io::read(path).unwrap_or_else(|e| panic!("read {what}: {e}"));
    assert_eq!(
        (img.nlines, img.nsamps, img.nbands),
        (250, 250, 4),
        "{what}: wrong shape"
    );
    assert_eq!(img.data, case1().data, "{what}: pixels differ from the ENVI original");

    let cfg = SegConfig {
        tols: vec![10.0],
        cm: 0.1,
        ..Default::default()
    }
    .with_n(&[15, 15, 100, 2500, 2500])
    .unwrap();

    let mut cap = Capture::default();
    run(img, &cfg, None, &mut cap).expect("segmentation");

    let (rp, rmap) = cap.rmap.unwrap();
    let (ap, armap) = cap.armap.unwrap();
    assert_eq!((rp, ap), (51, 58), "{what}: wrong pass counts");
    assert_eq!(
        rmap,
        std::fs::read(golden("test_3456/expected/proof/regmap.rmap.51")).unwrap(),
        "{what}: rmap does not match golden"
    );
    assert_eq!(
        armap,
        std::fs::read(golden("test_3456/expected/proof/regmap.armap.58")).unwrap(),
        "{what}: armap does not match golden"
    );
}

#[test]
fn tiff_roundtrip_reproduces_golden() {
    let img = case1();
    let path = tmp("case1.tif");
    {
        use tiff::encoder::{colortype, TiffEncoder};
        let f = std::fs::File::create(&path).unwrap();
        let mut enc = TiffEncoder::new(f).unwrap();
        // 4 bands -> a 4-sample TIFF. Photometric interpretation is irrelevant
        // here; the reader derives band count from bytes/pixel, which is what
        // multiband scientific imagery needs.
        enc.write_image::<colortype::CMYK8>(250, 250, img.data.as_u8().unwrap()).unwrap();
    }
    assert_reproduces_golden(&path, "TIFF");
}

#[test]
fn png_roundtrip_reproduces_golden() {
    let img = case1();
    let path = tmp("case1.png");
    {
        let f = std::fs::File::create(&path).unwrap();
        let w = std::io::BufWriter::new(f);
        let mut enc = png::Encoder::new(w, 250, 250);
        enc.set_color(png::ColorType::Rgba);
        enc.set_depth(png::BitDepth::Eight);
        let mut writer = enc.write_header().unwrap();
        writer.write_image_data(img.data.as_u8().unwrap()).unwrap();
    }
    assert_reproduces_golden(&path, "PNG");
}

/// Format detection must work off content, since the fixtures are inconsistent
/// about extensions (the real Case 1 input has none at all).
#[test]
fn detects_formats_by_content_not_extension() {
    use fast_segment::io::{detect, Format};
    let img = case1();

    let odd = tmp("no_extension_tiff");
    {
        use tiff::encoder::{colortype, TiffEncoder};
        let mut enc = TiffEncoder::new(std::fs::File::create(&odd).unwrap()).unwrap();
        enc.write_image::<colortype::CMYK8>(250, 250, img.data.as_u8().unwrap()).unwrap();
    }
    assert_eq!(detect(&odd).unwrap(), Format::Tiff);

    assert_eq!(detect(&golden("misc/temp_byte_bip")).unwrap(), Format::Envi);
    assert_eq!(
        detect(&golden("test_3456/input/test_3456.bip.ipw")).unwrap(),
        Format::Ipw
    );
}

/// 16-bit TIFF now reads, with the values intact rather than truncated to a
/// byte. This is the whole point of widening the input: a 1000 DN sample stays
/// 1000, not 232.
#[test]
fn reads_16bit_tiff() {
    let path = tmp("deep.tif");
    {
        use tiff::encoder::{colortype, TiffEncoder};
        let mut enc = TiffEncoder::new(std::fs::File::create(&path).unwrap()).unwrap();
        let data: Vec<u16> = (0..16 * 16u32).map(|i| 1000 + i as u16).collect();
        enc.write_image::<colortype::Gray16>(16, 16, &data).unwrap();
    }
    let img = fast_segment::io::read(&path).expect("16-bit TIFF should read");
    assert_eq!((img.nlines, img.nsamps, img.nbands), (16, 16, 1));
    let v = img.data.as_u16().expect("stored as u16");
    assert_eq!(v[0], 1000);
    assert_eq!(v[255], 1255);
}

/// 32-bit and floating-point samples are still refused: the algorithm's
/// distances and tolerances are integer DN.
#[test]
fn rejects_f32_tiff() {
    let path = tmp("float.tif");
    {
        use tiff::encoder::{colortype, TiffEncoder};
        let mut enc = TiffEncoder::new(std::fs::File::create(&path).unwrap()).unwrap();
        let data = vec![0.5f32; 16 * 16];
        enc.write_image::<colortype::Gray32Float>(16, 16, &data).unwrap();
    }
    let err = fast_segment::io::read(&path).expect_err("float TIFF should be rejected");
    assert!(
        err.to_string().contains("8- and 16-bit"),
        "unexpected: {err}"
    );
}
