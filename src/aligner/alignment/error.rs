use derive_more::From;
use std::ffi::{IntoStringError, NulError};
use std::fmt::{Display, Formatter};

#[derive(Debug, From)]
pub enum Error {
    NoStats(String),
    NoTable(String),
    NoStatsTable(String),
    NoRowCol(String),
    NoTrace(String),
    #[from]
    InvalidUTF8String(IntoStringError),
    #[from]
    InteriorNulByte(NulError),
    NoBandwidth,
}

impl Display for Error {
    fn fmt(&self, f: &mut Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{self:?}")
    }
}

impl std::error::Error for Error {}
