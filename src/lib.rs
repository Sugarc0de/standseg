//! A from-scratch Rust reimplementation of Harward & Woodcock's
//! nested-hierarchical region-growing image segmenter.
//!
//! See PLAN.md for the port strategy and the byte-exactness constraints.

pub mod config;
pub mod contig;
pub mod driver;
pub mod geo;
pub mod image;
pub mod io;
pub mod nbrset;
pub mod pixel;
pub mod region;
pub mod rng;
pub mod segment;
pub mod stage2;
pub mod vector;
