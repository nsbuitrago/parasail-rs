use derive_more::From;
use std::ffi::{CString, NulError};

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, From)]
pub enum Error {
    IncompatibleProfileConfig {
        vec_strategy: String,
    },
    InvalidSolutionWidth {
        solution_width: u8,
    },
    #[from]
    NulError(NulError),
    FunctionLookupFailed {
        fn_name: CString,
    },
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{self:?}")
    }
}

impl std::error::Error for Error {}
