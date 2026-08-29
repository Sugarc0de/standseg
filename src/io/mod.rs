//! Image input and output.
//!
//! Formats: ENVI (raw + `.hdr`) and IPW natively; TIFF and PNG to follow.
//! No libgdal -- the original's GDAL dependency covered exactly the two formats
//! we can write by hand. See PLAN.md section 9.

pub mod envi;
pub mod ipw;

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
        if ipw::sniff(&head[..n]) {
            return Ok(Format::Ipw);
        }
    }
    if let Some(ext) = path.extension().and_then(|s| s.to_str()) {
        match ext.to_lowercase().as_str() {
            "ipw" => return Ok(Format::Ipw),
            _ => {}
        }
    }
    // Anything with a readable ENVI sidecar is ENVI.
    let hdr = envi::header_path(path);
    if hdr.exists() {
        return Ok(Format::Envi);
    }
    Err(IoError::new(format!(
        "can't identify the format of {} -- not IPW, and no ENVI .hdr sidecar found \
         (looked for {})",
        path.display(),
        hdr.display()
    )))
}

pub fn read(path: &Path) -> Result<Image> {
    match detect(path)? {
        Format::Envi => envi::read(path),
        Format::Ipw => ipw::read(path),
    }
}
