//! The neighbour set used by `reg_nnbr`.
//!
//! **Insertion-ordered with linear dedup -- not a hash set.** `reg_nnbr` walks
//! this in insertion order, and that order decides which candidate first
//! establishes the running minimum, which decides how many `flip()` draws get
//! consumed. Reordering it would silently desync the RNG. PLAN.md 3.3.
//!
//! The C's `set.c` dedups by scanning backwards from the most recently added
//! entry; we do the same while the set is small. (Its 8-byte `SITEM` union and
//! the `case 4` branch that reads 8 bytes out of a 4-byte object -- flagged
//! `// OFFENDER` in the original -- are C genericity artifacts with no analogue
//! here.)
//!
//! The C also had a hard `MAX_NEIGHBORS` of 5000 and aborted the whole run on
//! the 5001st. That is gone: the list grows. What replaces it is a membership
//! side-table, because an unbounded backwards scan is quadratic. Past
//! `LINEAR_LIMIT` entries the dedup consults a hash set instead -- `items`
//! itself stays in exactly the same insertion order either way, so this cannot
//! move the RNG. See PLAN.md section 12.2.
//!
//! The side-table is boxed and allocated only on overflow. This struct is one
//! per scratch slot in the parallel sweep, hundreds of thousands of them, and
//! the sweep is the program's hot loop -- an inline `HashSet` would triple the
//! struct and cost about a tenth of total runtime for a branch that a normal
//! region never takes.

use std::collections::HashSet;

use crate::region::{RegionId, RegionList};

/// Below this many neighbours a backwards linear scan beats hashing. Typical
/// regions have a handful, so the hash path is close to unreachable in practice
/// -- it exists so a pathological region degrades instead of aborting.
pub const LINEAR_LIMIT: usize = 96;

#[derive(Default)]
pub struct NbrSet {
    /// Neighbour id and its squared centroid distance. The distance is filled
    /// in afterwards by `fill_dists`; the serial path ignores it and computes
    /// distances as it selects.
    items: Vec<(RegionId, f32)>,
    /// None until `items` reaches `LINEAR_LIMIT`, then authoritative.
    ///
    /// Boxed on purpose. `Segmenter::scratch` holds one `NbrSet` per region in
    /// a parallel batch, and almost none of them ever allocate the set, so the
    /// 8 bytes of a null pointer beat carrying a `HashSet` inline.
    #[allow(clippy::box_collection)]
    seen: Option<Box<HashSet<RegionId>>>,
}

impl NbrSet {
    pub fn new() -> Self {
        Self {
            items: Vec::with_capacity(64),
            seen: None,
        }
    }

    /// A non-allocating placeholder, so `reg_nnbr` can move the set out of the
    /// segmenter to satisfy the borrow checker without touching the allocator.
    pub const fn empty() -> Self {
        Self {
            items: Vec::new(),
            seen: None,
        }
    }

    #[inline]
    pub fn as_slice(&self) -> &[(RegionId, f32)] {
        &self.items
    }

    #[inline]
    pub fn clear(&mut self) {
        self.items.clear();
        // Drop the side-table rather than keep it: a run that once met a
        // 5000-neighbour region should not carry its hash set through the
        // hundreds of millions of ordinary ones that follow.
        self.seen = None;
    }

    /// Add `id` if it is not already present. Returns true if it was new.
    ///
    /// Unlike the C's `add_to_set`, this cannot fail.
    #[inline]
    pub fn add(&mut self, id: RegionId) -> bool {
        if self.items.len() < LINEAR_LIMIT {
            if self.items.iter().rev().any(|&(v, _)| v == id) {
                return false;
            }
            self.items.push((id, 0.0));
            if self.items.len() == LINEAR_LIMIT {
                self.start_hashing();
            }
            true
        } else {
            self.add_hashed(id)
        }
    }

    #[cold]
    #[inline(never)]
    fn start_hashing(&mut self) {
        let mut set = Box::new(HashSet::with_capacity(LINEAR_LIMIT * 2));
        set.extend(self.items.iter().map(|&(v, _)| v));
        self.seen = Some(set);
    }

    #[cold]
    #[inline(never)]
    fn add_hashed(&mut self, id: RegionId) -> bool {
        let seen = self
            .seen
            .as_mut()
            .expect("side-table exists past LINEAR_LIMIT");
        if seen.insert(id) {
            self.items.push((id, 0.0));
            true
        } else {
            false
        }
    }

    /// Fill each entry's squared centroid distance from `rid`.
    pub fn fill_dists(&mut self, rl: &RegionList, rid: RegionId) {
        for e in self.items.iter_mut() {
            e.1 = rl.dist2(rid, e.0);
        }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.items.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ids(s: &NbrSet) -> Vec<RegionId> {
        s.as_slice().iter().map(|&(v, _)| v).collect()
    }

    /// Insertion order is the load-bearing property; it must survive the
    /// switch from the linear scan to the hash set, and duplicates must still
    /// be rejected on both sides of the threshold.
    #[test]
    fn insertion_order_survives_the_hash_threshold() {
        let mut s = NbrSet::new();
        let n = LINEAR_LIMIT * 3;
        for i in 0..n {
            assert!(s.add(i as RegionId + 1), "id {i} should be new");
        }
        for i in 0..n {
            assert!(!s.add(i as RegionId + 1), "id {i} should be a duplicate");
        }
        assert_eq!(s.len(), n);
        assert_eq!(ids(&s), (1..=n as RegionId).collect::<Vec<_>>());
    }

    /// A duplicate straddling the threshold: added before the side-table
    /// exists, re-offered after it does.
    #[test]
    fn early_entries_are_still_deduped_after_the_switch() {
        let mut s = NbrSet::new();
        s.add(7);
        for i in 0..LINEAR_LIMIT + 10 {
            s.add(1000 + i as RegionId);
        }
        assert!(!s.add(7), "an early id must still be recognised as present");
        assert_eq!(s.as_slice()[0].0, 7);
    }

    #[test]
    fn clear_resets_both_paths() {
        let mut s = NbrSet::new();
        for i in 0..LINEAR_LIMIT + 5 {
            s.add(i as RegionId + 1);
        }
        s.clear();
        assert!(s.is_empty());
        assert!(s.add(1), "after clear, 1 is new again");
    }
}
