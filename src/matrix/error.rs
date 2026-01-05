use derive_more::From;
use std::{
    ffi::NulError,
    fmt::{Display, Formatter},
};

#[derive(Debug, From)]
pub enum Error {
    #[from]
    InteriorNulByte(NulError),
    FailedLookup(String),
    FileNotFound(String),
    NullMatrix,
    NotSquare,
    NotBuiltIn,
    InvalidIndex(i32, i32),
}

impl Display for Error {
    fn fmt(&self, f: &mut Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{self:?}")
    }
}

impl std::error::Error for Error {}
