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
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(rel)
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
            Phase::Stage2 => unreachable!("stage 2 was not asked for"),
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
    assert_eq!(
        img.data,
        case1().data,
        "{what}: pixels differ from the ENVI original"
    );

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
        enc.write_image::<colortype::CMYK8>(250, 250, img.data.as_u8().unwrap())
            .unwrap();
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
        enc.write_image::<colortype::CMYK8>(250, 250, img.data.as_u8().unwrap())
            .unwrap();
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

/// 32-bit float now *reads*, because stage-2 layers ship that way. It is still
/// refused by stage 1, and that half is pinned in `float_layer.rs` -- here we
/// only check the reader keeps the values as floats rather than rounding them
/// into an integer band, which is the failure that would be silent.
#[test]
fn reads_f32_tiff() {
    let path = tmp("float.tif");
    {
        use tiff::encoder::{colortype, TiffEncoder};
        let mut enc = TiffEncoder::new(std::fs::File::create(&path).unwrap()).unwrap();
        let data: Vec<f32> = (0..16 * 16u32).map(|i| 0.5 + i as f32 / 512.0).collect();
        enc.write_image::<colortype::Gray32Float>(16, 16, &data)
            .unwrap();
    }
    let img = fast_segment::io::read(&path).expect("float TIFF should read");
    assert_eq!((img.nlines, img.nsamps, img.nbands), (16, 16, 1));
    let v = img.data.as_f32().expect("stored as f32");
    assert_eq!(v[0], 0.5);
    assert_eq!(v[255], 0.5 + 255.0 / 512.0);
    assert!(img.data.is_float());
}

/// 64-bit float and the wider integers are still refused -- widening stopped at
/// what stage-2 layers actually are, and the message says where float belongs.
#[test]
fn rejects_f64_tiff() {
    let path = tmp("double.tif");
    {
        use tiff::encoder::{colortype, TiffEncoder};
        let mut enc = TiffEncoder::new(std::fs::File::create(&path).unwrap()).unwrap();
        let data = vec![0.5f64; 16 * 16];
        enc.write_image::<colortype::Gray64Float>(16, 16, &data)
            .unwrap();
    }
    let err = fast_segment::io::read(&path).expect_err("f64 TIFF should be rejected");
    assert!(
        err.to_string().contains("64-bit float") && err.to_string().contains("--stage2"),
        "unexpected: {err}"
    );
}

/// Cloud-optimized GeoTIFFs increasingly ship ZSTD-compressed, and the `tiff`
/// crate does *not* enable that codec by default -- it is a feature flag in
/// `Cargo.toml`. Dropping the flag would not break the build; it would make
/// real Landsat/Sentinel-derived COGs fail to open at runtime with
/// "compression method ZSTD is unsupported". This test is the guard on that
/// flag, so it writes its own file: the crate ships a ZSTD *decoder* but no
/// encoder, so we lay out a minimal single-strip TIFF by hand.
#[test]
fn reads_zstd_compressed_tiff() {
    const W: u32 = 16;
    const H: u32 = 16;
    let pixels: Vec<u16> = (0..W * H).map(|i| 1000 + i as u16).collect();
    let raw: Vec<u8> = pixels.iter().flat_map(|v| v.to_le_bytes()).collect();
    let strip = zstd::encode_all(&raw[..], 3).unwrap();

    // Tags must be written in ascending order; every value here fits the 4-byte
    // inline field, so there is no out-of-line value area to lay out.
    let tags: [(u16, u16, u32); 11] = [
        (256, 3, W),                  // ImageWidth
        (257, 3, H),                  // ImageLength
        (258, 3, 16),                 // BitsPerSample
        (259, 3, 0xC350),             // Compression = ZSTD
        (262, 3, 1),                  // Photometric = BlackIsZero
        (273, 4, 0),                  // StripOffsets -- patched below
        (277, 3, 1),                  // SamplesPerPixel
        (278, 3, H),                  // RowsPerStrip
        (279, 4, strip.len() as u32), // StripByteCounts
        (284, 3, 1),                  // PlanarConfiguration = chunky
        (339, 3, 1),                  // SampleFormat = unsigned integer
    ];
    let ifd_off = 8u32;
    let strip_off = ifd_off + 2 + 12 * tags.len() as u32 + 4;

    let mut f = Vec::new();
    f.extend_from_slice(b"II\x2a\x00");
    f.extend_from_slice(&ifd_off.to_le_bytes());
    f.extend_from_slice(&(tags.len() as u16).to_le_bytes());
    for (tag, kind, value) in tags {
        f.extend_from_slice(&tag.to_le_bytes());
        f.extend_from_slice(&kind.to_le_bytes());
        f.extend_from_slice(&1u32.to_le_bytes()); // count
        let v = if tag == 273 { strip_off } else { value };
        // A SHORT value sits in the low half of the field, not right-aligned.
        f.extend_from_slice(&v.to_le_bytes());
    }
    f.extend_from_slice(&0u32.to_le_bytes()); // no next IFD
    assert_eq!(f.len() as u32, strip_off);
    f.extend_from_slice(&strip);

    let path = tmp("zstd.tif");
    std::fs::write(&path, &f).unwrap();

    let img = fast_segment::io::read(&path).expect("ZSTD TIFF should read");
    assert_eq!((img.nlines, img.nsamps, img.nbands), (16, 16, 1));
    assert_eq!(img.data.as_u16().expect("stored as u16"), &pixels[..]);
}
