//! The neighbour set used by `reg_nnbr`.
//!
//! **Insertion-ordered with linear dedup -- not a hash set.** `reg_nnbr` walks
//! this in insertion order, and that order decides which candidate first
//! establishes the running minimum, which decides how many `flip()` draws get
//! consumed. Swapping in a HashSet would silently desync the RNG. PLAN.md 3.3.
//!
//! The C's `set.c` dedups by scanning backwards from the most recently added
//! entry; we do the same. (Its 8-byte `SITEM` union and the `case 4` branch that
//! reads 8 bytes out of a 4-byte object -- flagged `// OFFENDER` in the original
//! -- are C genericity artifacts with no analogue here.)

use crate::region::RegionId;

pub const MAX_NEIGHBORS: usize = 5000;

pub struct NbrSet {
    items: Vec<RegionId>,
    cap: usize,
}

impl NbrSet {
    pub fn new() -> Self {
        Self {
            items: Vec::with_capacity(64),
            cap: MAX_NEIGHBORS,
        }
    }

    /// A non-allocating placeholder, so `reg_nnbr` can move the set out of the
    /// segmenter to satisfy the borrow checker without touching the allocator.
    pub const fn empty() -> Self {
        Self {
            items: Vec::new(),
            cap: MAX_NEIGHBORS,
        }
    }

    #[inline]
    pub fn as_slice(&self) -> &[RegionId] {
        &self.items
    }

    #[inline]
    pub fn clear(&mut self) {
        self.items.clear();
    }

    /// Returns false only when the set is full, matching `add_to_set`.
    #[inline]
    pub fn add(&mut self, id: RegionId) -> bool {
        if self.items.iter().rev().any(|&v| v == id) {
            return true; // already present
        }
        if self.items.len() == self.cap {
            return false;
        }
        self.items.push(id);
        true
    }

    #[inline]
    pub fn iter(&self) -> impl Iterator<Item = RegionId> + '_ {
        self.items.iter().copied()
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.items.len()
    }
}

impl Default for NbrSet {
    fn default() -> Self {
        Self::new()
    }
}
