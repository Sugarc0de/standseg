//! Output files say what made them.
//!
//! IPW recorded `history = segment -t 10 -m .1 -n 15,15,100,2500,2500 ...` in
//! every image it wrote, and that record is the only reason the invocation
//! behind the golden fixtures was recoverable eleven years later. Our ENVI
//! output carried nothing, which was a regression against 1992.

use std::path::PathBuf;

use fast_segment::image::GeoRef;
use fast_segment::io::{self, Provenance};

fn tmp(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join("fast_segment_prov_test");
    std::fs::create_dir_all(&d).unwrap();
    d.join(name)
}

fn write_envi(path: &PathBuf, prov: &Provenance) {
    let rband: Vec<u32> = (0..16u32).collect();
    io::envi::write_region_map(path, &rband, 4, 4, 2, &GeoRef::default(), false, prov).unwrap();
}

#[test]
fn envi_header_records_the_command_and_stays_parseable() {
    let prov = Provenance::from_args(
        ["fast_segment", "-t", "10", "-m", ".1", "-o", "stands", "scene.tif"]
            .into_iter()
            .map(String::from),
    );
    let p = tmp("cmd.rmap.3");
    write_envi(&p, &prov);

    let hdr = std::fs::read_to_string(io::envi::header_path(&p)).unwrap();
    assert!(
        hdr.contains("history = {fast_segment -t 10 -m .1 -o stands scene.tif}"),
        "history line missing or mangled:\n{hdr}"
    );
    assert!(hdr.contains("software = {fast_segment "), "software line missing");

    // The added keys must not break the reader: a header we write is a header
    // we can read.
    let h = io::envi::read_header(&io::envi::header_path(&p)).expect("re-read our own header");
    assert_eq!((h.samples, h.lines, h.bands), (4, 4, 1));
    assert_eq!(h.data_type, 12);
}

/// A path with a space, a quote or a brace must not be able to corrupt the
/// header -- a `}` would close the block early and make the file unreadable.
#[test]
fn awkward_arguments_cannot_break_the_header() {
    let prov = Provenance::from_args(
        ["fast_segment", "-o", "my stands", "/data/{2014}/scene'1.tif"]
            .into_iter()
            .map(String::from),
    );
    let p = tmp("awkward.rmap.3");
    write_envi(&p, &prov);

    let hdr = std::fs::read_to_string(io::envi::header_path(&p)).unwrap();
    let line = hdr
        .lines()
        .find(|l| l.starts_with("history = "))
        .expect("no history line");
    assert!(line.ends_with('}'), "history line does not close: {line}");
    assert_eq!(line.matches('}').count(), 1, "stray brace in: {line}");
    assert!(line.contains("'my stands'"), "space not quoted: {line}");

    let h = io::envi::read_header(&io::envi::header_path(&p)).expect("header still parses");
    assert_eq!(h.bands, 1);
}

/// No timestamp anywhere: the same command twice gives byte-identical files.
/// Reproducibility is worth more here than knowing the hour of the run.
#[test]
fn output_is_deterministic() {
    let prov = Provenance::from_args(["fast_segment", "-t", "10"].into_iter().map(String::from));
    let (a, b) = (tmp("det_a.rmap.3"), tmp("det_b.rmap.3"));
    write_envi(&a, &prov);
    write_envi(&b, &prov);

    let ha = std::fs::read_to_string(io::envi::header_path(&a)).unwrap();
    let hb = std::fs::read_to_string(io::envi::header_path(&b)).unwrap();
    // The description and band name carry the file name, which differs by
    // design; everything else must match.
    let strip = |s: &str| {
        s.lines()
            .filter(|l| !l.contains("det_a") && !l.contains("det_b"))
            .collect::<Vec<_>>()
            .join("\n")
    };
    assert_eq!(strip(&ha), strip(&hb));
    assert_eq!(std::fs::read(&a).unwrap(), std::fs::read(&b).unwrap());
}

/// TIFF gets the same information in the tags the format already has for it.
#[test]
fn tiff_records_the_command_in_its_tags() {
    let prov = Provenance::from_args(
        ["fast_segment", "-t", "10", "-o", "stands"].into_iter().map(String::from),
    );
    let p = tmp("tagged.rmap.3.tif");
    let rband: Vec<u32> = (0..16u32).collect();
    io::tiff::write_region_map(&p, &rband, 4, 4, 2, &prov).unwrap();

    let raw = std::fs::read(&p).unwrap();
    let text = String::from_utf8_lossy(&raw);
    assert!(
        text.contains("fast_segment -t 10 -o stands"),
        "ImageDescription not written"
    );
    assert!(text.contains(env!("CARGO_PKG_VERSION")), "Software not written");

    // And it is still a readable TIFF.
    let img = io::read(&p).expect("re-read our own TIFF");
    assert_eq!((img.nlines, img.nsamps, img.nbands), (4, 4, 1));
}
