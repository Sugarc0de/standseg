//! In-memory image, band-interleaved-by-pixel, 8 bits per band.
//!
//! The original stores this as a dope vector (`uchar_t **`: an array of row
//! pointers into one contiguous block). That is a C ergonomics hack, not a space
//! saving -- it costs 8 bytes per row and an indirection per access. A flat
//! buffer with computed offsets is smaller and faster. See PLAN.md section 4, trick 1.

/// Georeferencing carried through from input to output where the format has it.
#[derive(Debug, Clone, Default)]
pub struct GeoRef {
    pub map_info: Option<String>,
    pub coord_sys: Option<String>,
    pub band_names: Vec<String>,
    pub description: Option<String>,
}

#[derive(Debug, Clone)]
pub struct Image {
    pub nlines: usize,
    pub nsamps: usize,
    pub nbands: usize,
    /// `nlines * nsamps * nbands` bytes, BIP order.
    pub data: Vec<u8>,
    pub geo: GeoRef,
}

impl Image {
    pub fn new(nlines: usize, nsamps: usize, nbands: usize) -> Self {
        Self {
            nlines,
            nsamps,
            nbands,
            data: vec![0u8; nlines * nsamps * nbands],
            geo: GeoRef::default(),
        }
    }

    #[inline]
    pub fn npixels(&self) -> usize {
        self.nlines * self.nsamps
    }

    /// The `nbands` values at (line, samp). Mirrors the C's `pixvec` macro.
    #[inline]
    pub fn pixel(&self, line: usize, samp: usize) -> &[u8] {
        let off = (line * self.nsamps + samp) * self.nbands;
        &self.data[off..off + self.nbands]
    }

    #[inline]
    pub fn pixel_at(&self, idx: usize) -> &[u8] {
        let off = idx * self.nbands;
        &self.data[off..off + self.nbands]
    }
}
