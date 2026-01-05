use derive_more::From;
use std::ffi::NulError;
use std::fmt::{Display, Formatter};

#[derive(Debug, From)]
pub enum Error {
    #[from]
    InteriorNulByte(NulError),
    NoBandwidth,
    #[from]
    Alignment(crate::aligner::alignment::Error),
}

impl Display for Error {
    fn fmt(&self, f: &mut Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{self:?}")
    }
}

impl std::error::Error for Error {}
