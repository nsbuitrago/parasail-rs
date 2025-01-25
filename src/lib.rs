mod aligner;
mod error;
mod matrix;
mod profile;

// Public API
pub use aligner::{Aligner, Gap};
pub use error::*;
pub use matrix::Matrix;
pub use profile::Profile;

#[derive(Debug, Clone, Copy)]
pub enum InstructionSet {
    SSE2,
    SSE41,
    AVX2,
    AltiVec,
    Neon,
}
