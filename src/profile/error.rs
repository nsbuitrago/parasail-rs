use derive_more::From;
use std::ffi::NulError;

use crate::InstructionSet;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, From)]
pub enum Error {
    InvalidSolutionWidth {
        solution_width: u8,
    },
    #[from]
    NullError(NulError),
    FunctionLookupFailed {
        use_stats: bool,
        instruction_set: Option<InstructionSet>,
        solution_width: Option<u8>,
    },
    NullProfile {
        msg: String,
    },
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{self:?}")
    }
}

impl std::error::Error for Error {}
