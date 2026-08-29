//! Image input and output.
//!
//! Formats: ENVI (raw + `.hdr`) and IPW natively; TIFF and PNG to follow.
//! No libgdal -- the original's GDAL dependency covered exactly the two formats
//! we can write by hand. See PLAN.md section 9.

pub mod envi;
pub mod ipw;
pub mod png;
pub mod tiff;

use std::fmt;
use std::path::Path;

use crate::image::Image;

#[derive(Debug)]
pub struct IoError(String);

impl IoError {
    pub fn new(msg: impl Into<String>) -> Self {
        IoError(msg.into())
    }
}

impl fmt::Display for IoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::error::Error for IoError {}

pub type Result<T> = std::result::Result<T, IoError>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Envi,
    Ipw,
    Tiff,
    Png,
}

/// Identify a file by content first, extension second.
///
/// Content wins because the fixtures are inconsistent about extensions: the
/// Case 1 ENVI input is named `temp_byte_bip` with no extension at all.
pub fn detect(path: &Path) -> Result<Format> {
    let mut head = [0u8; 16];
    {
        use std::io::Read;
        let mut f = std::fs::File::open(path)
            .map_err(|e| IoError::new(format!("can't open {}: {e}", path.display())))?;
        let n = f.read(&mut head).unwrap_or(0);
        let head = &head[..n];
        if ipw::sniff(head) {
            return Ok(Format::Ipw);
        }
        if png::sniff(head) {
            return Ok(Format::Png);
        }
        if tiff::sniff(head) {
            return Ok(Format::Tiff);
        }
    }
    if let Some(ext) = path.extension().and_then(|s| s.to_str()) {
        match ext.to_lowercase().as_str() {
            "ipw" => return Ok(Format::Ipw),
            "tif" | "tiff" => return Ok(Format::Tiff),
            "png" => return Ok(Format::Png),
            _ => {}
        }
    }
    // Anything with a readable ENVI sidecar is ENVI.
    let hdr = envi::header_path(path);
    if hdr.exists() {
        return Ok(Format::Envi);
    }
    Err(IoError::new(format!(
        "can't identify the format of {} -- not IPW, TIFF or PNG, and no ENVI .hdr \
         sidecar found (looked for {})",
        path.display(),
        hdr.display()
    )))
}

pub fn read(path: &Path) -> Result<Image> {
    Ok(read_with_nodata(path)?.0)
}

/// Read an image, along with any nodata value the format itself declares.
///
/// ENVI carries it as `data ignore value`, GeoTIFF as the `GDAL_NODATA` tag.
/// A `--nodata` on the command line overrides whatever comes back here.
pub fn read_with_nodata(path: &Path) -> Result<(Image, Option<f64>)> {
    match detect(path)? {
        Format::Envi => {
            let hdr = envi::read_header(&envi::header_path(path))?;
            let nd = hdr.data_ignore_value;
            Ok((envi::read(path)?, nd))
        }
        Format::Ipw => Ok((ipw::read(path)?, None)),
        Format::Tiff => {
            let r = tiff::read(path)?;
            Ok((r.image, r.nodata))
        }
        Format::Png => Ok((png::read(path)?, None)),
    }
}
