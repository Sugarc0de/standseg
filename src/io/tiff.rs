//! TIFF / GeoTIFF.
//!
//! Same sample widths as every other reader here: 8- and 16-bit integers. Bands
//! map to TIFF samples-per-pixel, so an RGB TIFF is a 3-band image and a 6-band
//! satellite stack is a 6-sample TIFF.

use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use tiff::decoder::{Decoder, DecodingResult};
use tiff::tags::Tag;

use crate::image::{GeoRef, Image, Samples};
use crate::io::{IoError, Result};

/// GDAL writes its nodata value into this private tag, as an ASCII string.
const GDAL_NODATA: u16 = 42113;

pub struct TiffRead {
    pub image: Image,
    /// Nodata declared by the file, if any.
    pub nodata: Option<f64>,
}

pub fn read(path: &Path) -> Result<TiffRead> {
    let f = File::open(path)
        .map_err(|e| IoError::new(format!("can't open {}: {e}", path.display())))?;
    let mut dec = Decoder::new(BufReader::new(f))
        .map_err(|e| IoError::new(format!("{}: not a readable TIFF: {e}", path.display())))?;

    let (w, h) = dec
        .dimensions()
        .map_err(|e| IoError::new(format!("{}: no dimensions: {e}", path.display())))?;

    let img = dec
        .read_image()
        .map_err(|e| IoError::new(format!("{}: can't decode: {e}", path.display())))?;

    let data = match img {
        DecodingResult::U8(v) => Samples::U8(v),
        DecodingResult::U16(v) => Samples::U16(v),
        DecodingResult::I16(v) => Samples::I16(v),
        other => {
            return Err(IoError::new(format!(
                "{}: TIFF samples are {}; this program segments 8- and 16-bit \
                 integer imagery only",
                path.display(),
                sample_kind(&other)
            )))
        }
    };

    let npix = w as usize * h as usize;
    if npix == 0 || data.len() % npix != 0 {
        return Err(IoError::new(format!(
            "{}: {} samples does not divide {}x{} pixels evenly",
            path.display(),
            data.len(),
            w,
            h
        )));
    }
    let nbands = data.len() / npix;

    // The tiff crate hands back interleaved samples, which is already BIP.
    let mut image = Image::from_samples(h as usize, w as usize, nbands, data);

    let nodata = dec
        .get_tag_ascii_string(Tag::Unknown(GDAL_NODATA))
        .ok()
        .and_then(|s| s.trim().parse::<f64>().ok());

    image.geo = GeoRef {
        description: path.file_name().map(|s| s.to_string_lossy().to_string()),
        ..Default::default()
    };

    Ok(TiffRead { image, nodata })
}

fn sample_kind(r: &DecodingResult) -> &'static str {
    match r {
        DecodingResult::U8(_) => "8-bit unsigned",
        DecodingResult::U16(_) => "16-bit unsigned",
        DecodingResult::U32(_) => "32-bit unsigned",
        DecodingResult::U64(_) => "64-bit unsigned",
        DecodingResult::I8(_) => "8-bit signed",
        DecodingResult::I16(_) => "16-bit signed",
        DecodingResult::I32(_) => "32-bit signed",
        DecodingResult::I64(_) => "64-bit signed",
        DecodingResult::F16(_) => "16-bit float",
        DecodingResult::F32(_) => "32-bit float",
        DecodingResult::F64(_) => "64-bit float",
    }
}

pub fn sniff(head: &[u8]) -> bool {
    head.starts_with(b"II\x2a\x00") || head.starts_with(b"MM\x00\x2a")
}

/// Write a single-band region map as a TIFF.
pub fn write_region_map(
    path: &Path,
    rband: &[u32],
    nlines: usize,
    nsamps: usize,
    nbytes: usize,
) -> Result<()> {
    use tiff::encoder::{colortype, TiffEncoder};

    let f = File::create(path)
        .map_err(|e| IoError::new(format!("can't create {}: {e}", path.display())))?;
    let mut enc = TiffEncoder::new(f)
        .map_err(|e| IoError::new(format!("{}: {e}", path.display())))?;

    let (w, h) = (nsamps as u32, nlines as u32);
    let n = nlines * nsamps;
    let res = match nbytes {
        1 => {
            let v: Vec<u8> = rband[..n].iter().map(|&x| x as u8).collect();
            enc.write_image::<colortype::Gray8>(w, h, &v)
        }
        2 => {
            let v: Vec<u16> = rband[..n].iter().map(|&x| x as u16).collect();
            enc.write_image::<colortype::Gray16>(w, h, &v)
        }
        4 => enc.write_image::<colortype::Gray32>(w, h, &rband[..n]),
        other => return Err(IoError::new(format!("bad region map pixel size {other}"))),
    };
    res.map_err(|e| IoError::new(format!("can't write {}: {e}", path.display())))
}
