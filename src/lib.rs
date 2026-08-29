//! A from-scratch Rust reimplementation of Harward & Woodcock's
//! nested-hierarchical region-growing image segmenter.
//!
//! See PLAN.md for the port strategy and the byte-exactness constraints.

pub mod image;
pub mod io;
pub mod rng;
