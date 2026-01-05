//! # Introduction
//!
//! This crate provides safe Rust bindings to
//! [Parasail](https://github.com/jeffdaily/parasail), a SIMD pairwise sequence
//! alignment library.
//!
//! For unsafe bindings, see the
//! [libparasail-sys](https://crates.io/crates/libparasail-sys) crate.
//!
//! # Usage
//!
//! ## Examples
//!
//! ### Basic usage
//! ```rust,no_run
//! use parasail_rs::{Aligner, Result};
//!
//! fn align() -> Result<()> {
//!     let query = b"ACGT";
//!     let reference = b"ACGT";
//!     let aligner = Aligner::new().build();
//!     let result = aligner.align(Some(query), reference)?;
//!     println!("Alignment Score: {}", result.get_score());
//!
//!     Ok(())
//! }
//! ```
//!
//! ### Using query profile
//!
//! When using striped or scan vectorization strategies, some performance may
//! be gained by reusing the query sequence. This can be done by creating a
//! query profile and reusing it for multiple alignments.
//!
//! ```rust,no_run
//! use parasail_rs::{Matrix, Aligner, Profile, Result};
//!
//! fn align() -> Result<()> {
//!
//!     let query = b"ACGT";
//!     let ref_1 = b"ACGTAACGTACA";
//!     let ref_2 = b"TGGCAAGGTAGA";
//!
//!     let use_stats = true;
//!     let query_profile = Profile::new(query, use_stats, &Matrix::default())?;
//!     let aligner = Aligner::new()
//!         .profile(query_profile)
//!         .build();
//!
//!     let result_1 = aligner.align(None, ref_1)?;
//!     let result_2 = aligner.align(None, ref_2)?;
//!
//!     println!("Score 1: {}", result_1.get_score());
//!     println!("Score 2: {}", result_2.get_score());
//!     Ok(())
//! }
//!

mod aligner;
mod error;
mod matrix;
mod profile;

pub use aligner::alignment::table::TraceFlags;
pub use aligner::Aligner;
pub use error::{Error, Result};
pub use matrix::Matrix;
pub use profile::Profile;

#[derive(Debug)]
pub enum SolutionWidth {
    Sat,
    Bit8,
    Bit16,
    Bit32,
    Bit64,
}

#[derive(Debug)]
pub enum InstructionSet {
    Best,
    SSE2,
    SSE41,
    AVX2,
    AltiVec,
    Neon,
}
