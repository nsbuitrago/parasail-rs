pub use crate::aligner::{Aligner, AlignerBuilder};
pub use crate::alignment::table::{Table, TraceFlags, TracebackTable};
pub use crate::alignment::{Alignment, SSWResult, Traceback};
pub use crate::error::{Error, Result};
pub use crate::matrix::Matrix;
pub use crate::profile::Profile;

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
