//! ENVI raw binary + `.hdr` sidecar.
//!
//! This is the format the golden fixtures use: `proof/regmap.armap.58` is a bare
//! 125000-byte raster with a separate text header. The original program wrote it
//! through GDAL's ENVI driver; we write the same bytes directly.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use crate::image::{GeoRef, Image};
use crate::io::{IoError, Result};

/// Parsed ENVI header. Only the fields the segmenter cares about.
#[derive(Debug, Clone, Default)]
pub struct EnviHeader {
    pub samples: usize,
    pub lines: usize,
    pub bands: usize,
    pub data_type: u32,
    pub interleave: String,
    pub byte_order: u32,
    pub header_offset: usize,
    pub data_ignore_value: Option<f64>,
    pub map_info: Option<String>,
    pub coord_sys: Option<String>,
    pub band_names: Vec<String>,
    pub description: Option<String>,
}

/// Split an ENVI header into `key = value` pairs, honouring `{ ... }` blocks
/// that may span lines.
fn parse_fields(text: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    let bytes: Vec<char> = text.chars().collect();
    let mut i = 0;
    // Skip the leading "ENVI" line.
    while i < bytes.len() && bytes[i] != '\n' {
        i += 1;
    }
    while i < bytes.len() {
        // Read a key up to '='.
        let start = i;
        while i < bytes.len() && bytes[i] != '=' && bytes[i] != '\n' {
            i += 1;
        }
        if i >= bytes.len() {
            break;
        }
        if bytes[i] == '\n' {
            i += 1;
            continue;
        }
        let key: String = bytes[start..i].iter().collect::<String>().trim().to_lowercase();
        i += 1; // past '='
        while i < bytes.len() && (bytes[i] == ' ' || bytes[i] == '\t') {
            i += 1;
        }
        let value: String = if i < bytes.len() && bytes[i] == '{' {
            i += 1;
            let vstart = i;
            let mut depth = 1;
            while i < bytes.len() && depth > 0 {
                match bytes[i] {
                    '{' => depth += 1,
                    '}' => depth -= 1,
                    _ => {}
                }
                if depth > 0 {
                    i += 1;
                }
            }
            let v: String = bytes[vstart..i].iter().collect();
            if i < bytes.len() {
                i += 1; // past '}'
            }
            v.trim().to_string()
        } else {
            let vstart = i;
            while i < bytes.len() && bytes[i] != '\n' {
                i += 1;
            }
            bytes[vstart..i].iter().collect::<String>().trim().to_string()
        };
        if !key.is_empty() {
            out.push((key, value));
        }
    }
    out
}

pub fn header_path(image_path: &Path) -> PathBuf {
    // ENVI convention allows both `img.hdr` and `img.bip.hdr`; the fixtures use
    // both (`test_3456.hdr` beside `test_3456.bip`, `regionmap.hdr` beside
    // `regionmap`). Try replacing the extension first, then appending.
    let with_replaced = image_path.with_extension("hdr");
    if with_replaced.exists() {
        return with_replaced;
    }
    let mut appended = image_path.as_os_str().to_owned();
    appended.push(".hdr");
    PathBuf::from(appended)
}

pub fn read_header(path: &Path) -> Result<EnviHeader> {
    let text = fs::read_to_string(path)
        .map_err(|e| IoError::new(format!("can't read ENVI header {}: {e}", path.display())))?;
    if !text.trim_start().starts_with("ENVI") {
        return Err(IoError::new(format!(
            "{} does not look like an ENVI header (no leading ENVI)",
            path.display()
        )));
    }
    let mut h = EnviHeader {
        interleave: "bsq".into(),
        ..Default::default()
    };
    for (k, v) in parse_fields(&text) {
        match k.as_str() {
            "samples" => h.samples = v.parse().unwrap_or(0),
            "lines" => h.lines = v.parse().unwrap_or(0),
            "bands" => h.bands = v.parse().unwrap_or(0),
            "data type" => h.data_type = v.parse().unwrap_or(0),
            "interleave" => h.interleave = v.to_lowercase(),
            "byte order" => h.byte_order = v.parse().unwrap_or(0),
            "header offset" => h.header_offset = v.parse().unwrap_or(0),
            "data ignore value" => h.data_ignore_value = v.parse().ok(),
            "map info" => h.map_info = Some(v),
            "coordinate system string" => h.coord_sys = Some(v),
            "description" => h.description = Some(v),
            "band names" => {
                h.band_names = v.split(',').map(|s| s.trim().to_string()).collect()
            }
            _ => {}
        }
    }
    if h.samples == 0 || h.lines == 0 || h.bands == 0 {
        return Err(IoError::new(format!(
            "{}: incomplete header (samples={}, lines={}, bands={})",
            path.display(),
            h.samples,
            h.lines,
            h.bands
        )));
    }
    Ok(h)
}

/// Read an ENVI image as uint8 BIP.
///
/// The original rejects anything but 8-bit data (`error("Image must be Byte
/// datatype")`), and so do we -- Case 2's int16 `_stack` is *supposed* to fail
/// here. See PLAN.md section 9.
pub fn read(path: &Path) -> Result<Image> {
    let hdr_path = header_path(path);
    let h = read_header(&hdr_path)?;

    if h.data_type != 1 {
        return Err(IoError::new(format!(
            "{}: ENVI data type {} is not 8-bit unsigned; this program segments \
             Byte imagery only (as did the original)",
            path.display(),
            h.data_type
        )));
    }

    let raw = fs::read(path)
        .map_err(|e| IoError::new(format!("can't read {}: {e}", path.display())))?;
    let want = h.lines * h.samples * h.bands;
    let avail = raw.len().saturating_sub(h.header_offset);
    if avail < want {
        return Err(IoError::new(format!(
            "{}: short read -- header says {} bytes ({}x{}x{}), file has {}",
            path.display(),
            want,
            h.samples,
            h.lines,
            h.bands,
            avail
        )));
    }
    let raw = &raw[h.header_offset..h.header_offset + want];

    let mut img = Image::new(h.lines, h.samples, h.bands);
    let (nl, ns, nb) = (h.lines, h.samples, h.bands);
    match h.interleave.as_str() {
        "bip" => img.data.copy_from_slice(raw),
        "bsq" => {
            for b in 0..nb {
                let base = b * nl * ns;
                for p in 0..nl * ns {
                    img.data[p * nb + b] = raw[base + p];
                }
            }
        }
        "bil" => {
            for l in 0..nl {
                for b in 0..nb {
                    let base = (l * nb + b) * ns;
                    for s in 0..ns {
                        img.data[(l * ns + s) * nb + b] = raw[base + s];
                    }
                }
            }
        }
        other => {
            return Err(IoError::new(format!(
                "{}: unsupported ENVI interleave '{other}'",
                path.display()
            )))
        }
    }

    img.geo = GeoRef {
        map_info: h.map_info,
        coord_sys: h.coord_sys,
        band_names: h.band_names,
        description: h.description,
    };
    Ok(img)
}

/// Write a single-band region map: raw little-endian pixels plus a `.hdr`.
///
/// `nbytes` (1, 2 or 4) is chosen by the caller exactly as `GDAL_write_image`
/// does, from the bit width of the largest region id.
pub fn write_region_map(
    path: &Path,
    rband: &[u32],
    nlines: usize,
    nsamps: usize,
    nbytes: usize,
    geo: &GeoRef,
    masked_present: bool,
) -> Result<()> {
    let mut out = Vec::with_capacity(nlines * nsamps * nbytes);
    for &v in rband.iter().take(nlines * nsamps) {
        let le = v.to_le_bytes();
        out.extend_from_slice(&le[..nbytes]);
    }
    fs::write(path, &out)
        .map_err(|e| IoError::new(format!("can't write {}: {e}", path.display())))?;

    let data_type = match nbytes {
        1 => 1,  // uint8
        2 => 12, // uint16
        4 => 13, // uint32
        n => return Err(IoError::new(format!("bad region map pixel size {n}"))),
    };

    let name = path
        .file_name()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_default();
    let mut hdr = String::from("ENVI\n");
    hdr.push_str(&format!("description = {{\n{name}}}\n"));
    hdr.push_str(&format!("samples = {nsamps}\n"));
    hdr.push_str(&format!("lines   = {nlines}\n"));
    hdr.push_str("bands   = 1\n");
    hdr.push_str("header offset = 0\n");
    hdr.push_str("file type = ENVI Standard\n");
    hdr.push_str(&format!("data type = {data_type}\n"));
    hdr.push_str("interleave = bsq\n");
    hdr.push_str("byte order = 0\n");
    if masked_present {
        // Region 0 is the artificial region holding masked / nodata pixels, so
        // water and non-treed area round-trip as nodata rather than as a stand.
        hdr.push_str("data ignore value = 0\n");
    }
    if let Some(mi) = &geo.map_info {
        hdr.push_str(&format!("map info = {{{mi}}}\n"));
    }
    if let Some(cs) = &geo.coord_sys {
        hdr.push_str(&format!("coordinate system string = {{{cs}}}\n"));
    }
    hdr.push_str("band names = {\nBand 1}\n");

    let hp = header_path(path);
    fs::write(&hp, hdr)
        .map_err(|e| IoError::new(format!("can't write {}: {e}", hp.display())))?;
    Ok(())
}
