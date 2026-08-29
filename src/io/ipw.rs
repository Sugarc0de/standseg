//! IPW (Image Processing Workbench) images.
//!
//! Header is a sequence of plain-text records:
//!
//! ```text
//! !<header> basic_image_i -1 $Revision: 1.10 $
//! byteorder = 01234567
//! nlines = 250
//! nsamps = 250
//! nbands = 4
//! !<header> basic_image 0 $Revision: 1.10 $
//! bytes = 1
//! bits = 8
//! ...
//! !<header> image -1 $Revision: 1.5 $
//! <form feed><newline><raw pixels>
//! ```
//!
//! Pixels are **byte-aligned** at `bytes` per pixel and masked to `bits` -- not
//! bit-packed (see `libipw/pixio/pvwrite.c`, which copies whole bytes). Bands
//! are interleaved by pixel.

use std::fs;
use std::path::Path;

use crate::image::{GeoRef, Image, Samples};
use crate::io::{IoError, Result};

#[derive(Debug, Clone, Default)]
pub struct IpwHeader {
    pub nlines: usize,
    pub nsamps: usize,
    pub nbands: usize,
    /// Per-band bytes per pixel.
    pub bytes: Vec<usize>,
    /// Per-band significant bits.
    pub bits: Vec<usize>,
    pub byteorder: Option<String>,
    pub history: Vec<String>,
    /// Offset of the first pixel byte.
    pub data_offset: usize,
}

/// Locate the end of the header: the form feed that closes the last record.
fn find_data_offset(raw: &[u8]) -> Option<usize> {
    let ff = raw.iter().position(|&b| b == 0x0c)?;
    // A newline conventionally follows the form feed; skip it if present.
    if raw.get(ff + 1) == Some(&b'\n') {
        Some(ff + 2)
    } else {
        Some(ff + 1)
    }
}

pub fn parse_header(raw: &[u8]) -> Result<IpwHeader> {
    let data_offset = find_data_offset(raw)
        .ok_or_else(|| IoError::new("not an IPW file: no form feed terminating the header"))?;
    let text = String::from_utf8_lossy(&raw[..data_offset]);

    let mut h = IpwHeader {
        data_offset,
        ..Default::default()
    };
    // Which band the current `basic_image N` record describes; -1 is the
    // image-wide `basic_image_i` record.
    let mut cur_band: i32 = -1;

    for line in text.lines() {
        let line = line.trim_end_matches('\u{c}');
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(rest) = line.strip_prefix("!<header>") {
            let mut it = rest.split_whitespace();
            let name = it.next().unwrap_or("");
            let band: i32 = it.next().and_then(|s| s.parse().ok()).unwrap_or(-1);
            cur_band = if name == "basic_image" { band } else { -1 };
            continue;
        }
        let Some((k, v)) = line.split_once('=') else {
            continue;
        };
        let (k, v) = (k.trim(), v.trim());
        match k {
            "nlines" => h.nlines = v.parse().unwrap_or(0),
            "nsamps" => h.nsamps = v.parse().unwrap_or(0),
            "nbands" => h.nbands = v.parse().unwrap_or(0),
            "byteorder" => h.byteorder = Some(v.to_string()),
            "history" => h.history.push(v.to_string()),
            "bytes" | "bits" => {
                let n: usize = v.parse().unwrap_or(0);
                if cur_band >= 0 {
                    let b = cur_band as usize;
                    let target = if k == "bytes" { &mut h.bytes } else { &mut h.bits };
                    if target.len() <= b {
                        target.resize(b + 1, 0);
                    }
                    target[b] = n;
                }
            }
            _ => {}
        }
    }

    if h.nlines == 0 || h.nsamps == 0 || h.nbands == 0 {
        return Err(IoError::new(format!(
            "IPW header incomplete: nlines={}, nsamps={}, nbands={}",
            h.nlines, h.nsamps, h.nbands
        )));
    }
    if h.bytes.len() < h.nbands || h.bits.len() < h.nbands {
        return Err(IoError::new(format!(
            "IPW header declares {} bands but only {} bytes/{} bits records",
            h.nbands,
            h.bytes.len(),
            h.bits.len()
        )));
    }
    Ok(h)
}

/// Read an IPW image as BIP.
///
/// 1- and 2-byte bands are supported; all bands must be the same width, since
/// the in-memory image is one sample type throughout. Values are masked to the
/// band's declared `bits`, as `libipw` does.
pub fn read(path: &Path) -> Result<Image> {
    let raw = fs::read(path)
        .map_err(|e| IoError::new(format!("can't read {}: {e}", path.display())))?;
    let h = parse_header(&raw)?;

    let nbytes = h.bytes[0];
    for b in 0..h.nbands {
        if h.bytes[b] != nbytes {
            return Err(IoError::new(format!(
                "{}: band 0 is {} bytes/pixel but band {} is {}; mixed sample widths \
                 are not supported",
                path.display(),
                nbytes,
                b,
                h.bytes[b]
            )));
        }
        if h.bits[b] == 0 || h.bits[b] > 8 * h.bytes[b] {
            return Err(IoError::new(format!(
                "{}: band {} declares {} bits in {} bytes, which is impossible",
                path.display(),
                b,
                h.bits[b],
                h.bytes[b]
            )));
        }
    }
    if nbytes != 1 && nbytes != 2 {
        return Err(IoError::new(format!(
            "{}: {} bytes per pixel; this program segments 8- and 16-bit imagery only",
            path.display(),
            nbytes
        )));
    }

    let npix = h.nlines * h.nsamps;
    let want = npix * h.nbands * nbytes;
    let avail = raw.len() - h.data_offset;
    if avail < want {
        return Err(IoError::new(format!(
            "{}: short read -- header says {want} pixel bytes, file has {avail}",
            path.display()
        )));
    }
    let pixels = &raw[h.data_offset..h.data_offset + want];

    // IPW records the writing machine's byte order as a permutation string;
    // "01234567" is MSB-first. Anything else means the bytes were written in a
    // little-endian order.
    let be = h
        .byteorder
        .as_deref()
        .map_or(true, |s| s.trim().starts_with('0'));

    let data = if nbytes == 1 {
        let mut v = pixels.to_vec();
        // Mask each band to its declared bits.
        if h.bits[..h.nbands].iter().any(|&b| b < 8) {
            let masks: Vec<u8> = (0..h.nbands)
                .map(|b| (((1u16 << h.bits[b]) - 1) & 0xff) as u8)
                .collect();
            for (i, x) in v.iter_mut().enumerate() {
                *x &= masks[i % h.nbands];
            }
        }
        Samples::U8(v)
    } else {
        let masks: Vec<u16> = (0..h.nbands)
            .map(|b| {
                if h.bits[b] >= 16 {
                    u16::MAX
                } else {
                    (1u16 << h.bits[b]) - 1
                }
            })
            .collect();
        let v: Vec<u16> = pixels
            .chunks_exact(2)
            .enumerate()
            .map(|(i, c)| {
                let w = [c[0], c[1]];
                let x = if be {
                    u16::from_be_bytes(w)
                } else {
                    u16::from_le_bytes(w)
                };
                x & masks[i % h.nbands]
            })
            .collect();
        Samples::U16(v)
    };

    let mut img = Image::from_samples(h.nlines, h.nsamps, h.nbands, data);
    img.geo = GeoRef {
        description: h.history.first().cloned(),
        ..Default::default()
    };
    Ok(img)
}

/// True if the file starts with an IPW header record.
pub fn sniff(raw: &[u8]) -> bool {
    raw.starts_with(b"!<header>")
}
