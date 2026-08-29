//! In-memory image, band-interleaved-by-pixel.
//!
//! The original stores this as a dope vector (`uchar_t **`: an array of row
//! pointers into one contiguous block). That is a C ergonomics hack, not a space
//! saving -- it costs 8 bytes per row and an indirection per access. A flat
//! buffer with computed offsets is smaller and faster. See PLAN.md section 4, trick 1.
//!
//! # Sample width
//!
//! The 1992 program was uint8-only, because that is what an 8-bit TM scene was.
//! Landsat 8/9 and Sentinel-2 are 12-bit data delivered in 16-bit containers, so
//! forcing them down to 8 bits throws away radiometry *before* segmentation and
//! changes the answer. We therefore carry three sample widths -- `u8`, `u16` and
//! `i16` (the container Landsat Collection 2 surface reflectance ships in).
//!
//! Nothing downstream of Phase 0 sees a pixel: the region centroids are f32 and
//! the image is freed as soon as the region list exists. So widening the input
//! is confined to this file, `pixel.rs` and `RegionList::from_pixel`, and the
//! `u8` path keeps its exact arithmetic -- which the golden fixtures pin.

/// A pixel sample type the segmenter can read.
pub trait Sample: Copy + PartialEq + Send + Sync + 'static {
    /// Name used in diagnostics.
    const KIND: &'static str;
    const MIN_VALUE: i64;
    const MAX_VALUE: i64;
    fn to_i64(self) -> i64;
    fn to_f32(self) -> f32;
}

impl Sample for u8 {
    const KIND: &'static str = "8-bit unsigned";
    const MIN_VALUE: i64 = 0;
    const MAX_VALUE: i64 = u8::MAX as i64;
    #[inline]
    fn to_i64(self) -> i64 {
        self as i64
    }
    #[inline]
    fn to_f32(self) -> f32 {
        self as f32
    }
}

impl Sample for u16 {
    const KIND: &'static str = "16-bit unsigned";
    const MIN_VALUE: i64 = 0;
    const MAX_VALUE: i64 = u16::MAX as i64;
    #[inline]
    fn to_i64(self) -> i64 {
        self as i64
    }
    #[inline]
    fn to_f32(self) -> f32 {
        self as f32
    }
}

impl Sample for i16 {
    const KIND: &'static str = "16-bit signed";
    const MIN_VALUE: i64 = i16::MIN as i64;
    const MAX_VALUE: i64 = i16::MAX as i64;
    #[inline]
    fn to_i64(self) -> i64 {
        self as i64
    }
    #[inline]
    fn to_f32(self) -> f32 {
        self as f32
    }
}

/// Pixel storage, one variant per supported sample width. BIP order throughout:
/// `nlines * nsamps * nbands` samples, bands adjacent within a pixel.
#[derive(Debug, Clone, PartialEq)]
pub enum Samples {
    U8(Vec<u8>),
    U16(Vec<u16>),
    I16(Vec<i16>),
}

impl Samples {
    pub fn len(&self) -> usize {
        match self {
            Samples::U8(v) => v.len(),
            Samples::U16(v) => v.len(),
            Samples::I16(v) => v.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Bytes of memory per sample -- what the memory report needs.
    pub fn bytes_per_sample(&self) -> usize {
        match self {
            Samples::U8(_) => 1,
            Samples::U16(_) | Samples::I16(_) => 2,
        }
    }

    pub fn kind(&self) -> &'static str {
        match self {
            Samples::U8(_) => <u8 as Sample>::KIND,
            Samples::U16(_) => <u16 as Sample>::KIND,
            Samples::I16(_) => <i16 as Sample>::KIND,
        }
    }

    pub fn as_u8(&self) -> Option<&[u8]> {
        match self {
            Samples::U8(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_u16(&self) -> Option<&[u16]> {
        match self {
            Samples::U16(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_i16(&self) -> Option<&[i16]> {
        match self {
            Samples::I16(v) => Some(v),
            _ => None,
        }
    }

    /// Inclusive range a value must sit in to be a valid sample of this type.
    pub fn value_range(&self) -> (i64, i64) {
        match self {
            Samples::U8(_) => (<u8 as Sample>::MIN_VALUE, <u8 as Sample>::MAX_VALUE),
            Samples::U16(_) => (<u16 as Sample>::MIN_VALUE, <u16 as Sample>::MAX_VALUE),
            Samples::I16(_) => (<i16 as Sample>::MIN_VALUE, <i16 as Sample>::MAX_VALUE),
        }
    }
}

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
    pub data: Samples,
    pub geo: GeoRef,
}

/// A typed, borrowed view of the pixel buffer. Phase 0 is generic over this so
/// each sample width gets its own specialised code, rather than paying a match
/// per pixel access.
pub struct Raster<'a, T> {
    pub nlines: usize,
    pub nsamps: usize,
    pub nbands: usize,
    pub data: &'a [T],
}

impl<'a, T: Sample> Raster<'a, T> {
    /// The `nbands` values at (line, samp). Mirrors the C's `pixvec` macro.
    #[inline]
    pub fn pixel(&self, line: usize, samp: usize) -> &[T] {
        let off = (line * self.nsamps + samp) * self.nbands;
        &self.data[off..off + self.nbands]
    }

    #[inline]
    pub fn pixel_at(&self, idx: usize) -> &[T] {
        let off = idx * self.nbands;
        &self.data[off..off + self.nbands]
    }

    #[inline]
    pub fn npixels(&self) -> usize {
        self.nlines * self.nsamps
    }
}

impl Image {
    /// An 8-bit image of zeros. 8-bit stays the default because it is what the
    /// fixtures are and what most callers construct.
    pub fn new(nlines: usize, nsamps: usize, nbands: usize) -> Self {
        Self::from_samples(
            nlines,
            nsamps,
            nbands,
            Samples::U8(vec![0u8; nlines * nsamps * nbands]),
        )
    }

    pub fn from_samples(nlines: usize, nsamps: usize, nbands: usize, data: Samples) -> Self {
        debug_assert_eq!(data.len(), nlines * nsamps * nbands);
        Self {
            nlines,
            nsamps,
            nbands,
            data,
            geo: GeoRef::default(),
        }
    }

    #[inline]
    pub fn npixels(&self) -> usize {
        self.nlines * self.nsamps
    }

    /// Bytes the pixel buffer occupies.
    pub fn nbytes(&self) -> usize {
        self.data.len() * self.data.bytes_per_sample()
    }

    /// Run `f` against whichever typed view the image actually holds.
    ///
    /// This is the single dispatch point: everything below it is monomorphic.
    pub fn with_raster<R>(&self, f: impl FnOnce(RasterRef<'_>) -> R) -> R {
        let (nl, ns, nb) = (self.nlines, self.nsamps, self.nbands);
        match &self.data {
            Samples::U8(v) => f(RasterRef::U8(Raster {
                nlines: nl,
                nsamps: ns,
                nbands: nb,
                data: v,
            })),
            Samples::U16(v) => f(RasterRef::U16(Raster {
                nlines: nl,
                nsamps: ns,
                nbands: nb,
                data: v,
            })),
            Samples::I16(v) => f(RasterRef::I16(Raster {
                nlines: nl,
                nsamps: ns,
                nbands: nb,
                data: v,
            })),
        }
    }

    /// Every pixel whose bands match `nd` is marked 0 in `mask`.
    ///
    /// `any` picks the "nodata if any band matches" reading; the default is that
    /// all bands must match. A value outside the sample type's range simply
    /// matches nothing.
    pub fn apply_nodata(&self, nd: i64, any: bool, mask: &mut [u8]) {
        fn scan<T: Sample>(r: &Raster<'_, T>, nd: i64, any: bool, mask: &mut [u8]) {
            for p in 0..r.npixels() {
                let px = r.pixel_at(p);
                let hit = if any {
                    px.iter().any(|s| s.to_i64() == nd)
                } else {
                    px.iter().all(|s| s.to_i64() == nd)
                };
                if hit {
                    mask[p] = 0;
                }
            }
        }
        self.with_raster(|r| match r {
            RasterRef::U8(r) => scan(&r, nd, any, mask),
            RasterRef::U16(r) => scan(&r, nd, any, mask),
            RasterRef::I16(r) => scan(&r, nd, any, mask),
        })
    }

    /// Flatten a single-band image to a 0/1 mask: any nonzero sample is valid.
    pub fn to_mask(&self) -> Vec<u8> {
        fn flat<T: Sample>(r: &Raster<'_, T>) -> Vec<u8> {
            r.data
                .iter()
                .map(|s| u8::from(s.to_i64() != 0))
                .collect()
        }
        self.with_raster(|r| match r {
            RasterRef::U8(r) => flat(&r),
            RasterRef::U16(r) => flat(&r),
            RasterRef::I16(r) => flat(&r),
        })
    }
}

/// The three typed views `Image::with_raster` can hand back.
pub enum RasterRef<'a> {
    U8(Raster<'a, u8>),
    U16(Raster<'a, u16>),
    I16(Raster<'a, i16>),
}
