//! PlanarConfiguration = 2.
//!
//! A multi-band TIFF can store its samples two ways: interleaved per pixel
//! (`PlanarConfiguration = 1`, "chunky") or one whole plane after another
//! (`= 2`, "planar"). GDAL writes both, and the NTEMS species probability
//! stacks are planar.
//!
//! The `tiff` crate's `read_image` decodes only the *first plane*, and a 6-band
//! planar file then divides evenly into what looks like a perfectly good 1-band
//! image. Nothing errors; the segmentation just silently runs on band 1. This
//! pins the fix: both layouts must produce the same image.

use fast_segment::image::Samples;

const W: usize = 4;
const H: usize = 3;
const SPP: usize = 2;

/// Hand-built so the test owns the bytes: the crate's encoder writes chunky
/// only, and the file this is really about lives on an external drive.
fn tiff(planar: bool) -> Vec<u8> {
    // Band 0 counts up from 10, band 1 from 100, so a plane mix-up is obvious.
    let plane0: Vec<u8> = (0..W * H).map(|i| 10 + i as u8).collect();
    let plane1: Vec<u8> = (0..W * H).map(|i| 100 + i as u8).collect();

    const IFD: u32 = 8;
    const NENT: u16 = 12;
    let data_at = IFD + 2 + u32::from(NENT) * 12 + 4; // 134

    let mut e: Vec<(u16, u16, u32, [u8; 4])> = Vec::new(); // tag, type, count, value
    let short1 = |v: u16| {
        let b = v.to_le_bytes();
        [b[0], b[1], 0, 0]
    };
    let short2 = |a: u16, b: u16| {
        let (a, b) = (a.to_le_bytes(), b.to_le_bytes());
        [a[0], a[1], b[0], b[1]]
    };
    e.push((256, 3, 1, short1(W as u16))); // ImageWidth
    e.push((257, 3, 1, short1(H as u16))); // ImageLength
    e.push((258, 3, 2, short2(8, 8))); // BitsPerSample
    e.push((259, 3, 1, short1(1))); // Compression: none
    e.push((262, 3, 1, short1(1))); // Photometric: BlackIsZero
    if planar {
        // One strip per plane, back to back.
        e.push((273, 3, 2, short2(data_at as u16, data_at as u16 + 12)));
    } else {
        e.push((273, 3, 1, short1(data_at as u16)));
    }
    e.push((277, 3, 1, short1(SPP as u16))); // SamplesPerPixel
    e.push((278, 3, 1, short1(H as u16))); // RowsPerStrip: the whole image
    if planar {
        e.push((279, 3, 2, short2(12, 12))); // StripByteCounts
    } else {
        e.push((279, 3, 1, short1(24)));
    }
    e.push((284, 3, 1, short1(if planar { 2 } else { 1 }))); // PlanarConfiguration
    e.push((338, 3, 1, short1(0))); // ExtraSamples: unspecified
    e.push((339, 3, 2, short2(1, 1))); // SampleFormat: unsigned integer
    assert_eq!(e.len(), NENT as usize);

    let mut f = Vec::new();
    f.extend_from_slice(b"II");
    f.extend_from_slice(&42u16.to_le_bytes());
    f.extend_from_slice(&IFD.to_le_bytes());
    f.extend_from_slice(&NENT.to_le_bytes());
    for (tag, ty, count, val) in &e {
        f.extend_from_slice(&tag.to_le_bytes());
        f.extend_from_slice(&ty.to_le_bytes());
        f.extend_from_slice(&count.to_le_bytes());
        f.extend_from_slice(val);
    }
    f.extend_from_slice(&0u32.to_le_bytes()); // no next IFD
    assert_eq!(f.len(), data_at as usize);

    if planar {
        f.extend_from_slice(&plane0);
        f.extend_from_slice(&plane1);
    } else {
        for i in 0..W * H {
            f.push(plane0[i]);
            f.push(plane1[i]);
        }
    }
    f
}

fn read(planar: bool) -> fast_segment::image::Image {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("build/out/tiff_planar");
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join(if planar { "planar.tif" } else { "chunky.tif" });
    std::fs::write(&p, tiff(planar)).unwrap();
    fast_segment::io::read(&p).unwrap_or_else(|e| panic!("read {}: {e}", p.display()))
}

/// The two layouts describe the same image, so they must read as the same
/// image. Before the fix the planar one read as 1 band of 12 pixels.
#[test]
fn planar_and_chunky_read_identically() {
    let c = read(false);
    let p = read(true);
    assert_eq!(
        (p.nlines, p.nsamps, p.nbands),
        (H, W, SPP),
        "planar TIFF read as {} bands -- only the first plane was decoded",
        p.nbands
    );
    assert_eq!((c.nlines, c.nsamps, c.nbands), (H, W, SPP));
    assert_eq!(
        c.data, p.data,
        "planar samples were not deinterleaved to BIP"
    );

    // And the values are actually right, not merely equal to each other.
    let want: Vec<u8> = (0..W * H)
        .flat_map(|i| [10 + i as u8, 100 + i as u8])
        .collect();
    assert_eq!(p.data, Samples::U8(want), "wrong sample order");
}

/// A band that exists must reach the segmenter. This is the property the bug
/// actually violated: `-B 5` on a 6-band planar stack addressed a band the
/// program had thrown away.
#[test]
fn every_band_survives_the_read() {
    let p = read(true);
    let s = p.data.as_u8().unwrap();
    for b in 0..SPP {
        let band: Vec<u8> = (0..W * H).map(|i| s[i * SPP + b]).collect();
        let base = if b == 0 { 10u8 } else { 100u8 };
        assert_eq!(
            band,
            (0..W * H).map(|i| base + i as u8).collect::<Vec<_>>(),
            "band {b} is wrong"
        );
    }
}
