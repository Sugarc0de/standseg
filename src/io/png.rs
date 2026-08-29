//! PNG.
//!
//! 8-bit only. Channel count becomes the band count: grayscale is 1 band, RGB
//! is 3, RGBA is 4 -- note that an alpha channel is read as an ordinary band and
//! will take part in the spectral distance, which is almost never what you want.
//! Use `--nodata` or a mask for transparency instead.

use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use crate::image::{GeoRef, Image};
use crate::io::{IoError, Result};

pub fn read(path: &Path) -> Result<Image> {
    let f = File::open(path)
        .map_err(|e| IoError::new(format!("can't open {}: {e}", path.display())))?;
    let dec = png::Decoder::new(BufReader::new(f));
    let mut reader = dec
        .read_info()
        .map_err(|e| IoError::new(format!("{}: not a readable PNG: {e}", path.display())))?;

    let info = reader.info();
    if info.bit_depth != png::BitDepth::Eight {
        return Err(IoError::new(format!(
            "{}: PNG is {:?} bits per sample; this program segments 8-bit \
             imagery only (as did the original)",
            path.display(),
            info.bit_depth
        )));
    }
    let (w, h) = (info.width as usize, info.height as usize);

    let mut buf = vec![0u8; reader.output_buffer_size().unwrap_or(0)];
    let out = reader
        .next_frame(&mut buf)
        .map_err(|e| IoError::new(format!("{}: can't decode: {e}", path.display())))?;
    let nbands = out.color_type.samples();
    buf.truncate(out.buffer_size());

    if buf.len() != w * h * nbands {
        return Err(IoError::new(format!(
            "{}: decoded {} bytes, expected {}x{}x{}",
            path.display(),
            buf.len(),
            w,
            h,
            nbands
        )));
    }

    // PNG rows are already interleaved by pixel, i.e. BIP.
    let mut image = Image::new(h, w, nbands);
    image.data = buf;
    image.geo = GeoRef {
        description: path.file_name().map(|s| s.to_string_lossy().to_string()),
        ..Default::default()
    };
    Ok(image)
}

pub fn sniff(head: &[u8]) -> bool {
    head.starts_with(b"\x89PNG\r\n\x1a\n")
}
