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
    EmptyAlphabet,
    EmptyMatrixName,
    InvalidScores {
        match_score: i32,
        mismatch_score: i32,
    },
    NullMatrix,
    NotSquare,
    NotBuiltIn,
    InvalidIndex(i32, i32),
    InvalidPSSMValues {
        alphabet_len: usize,
        rows: i32,
        expected: usize,
        actual: usize,
    },
    EmptyPSSMAlphabet,
    InvalidPSSMRows(i32),
    PSSMTooLarge {
        alphabet_len: usize,
        rows: i32,
    },
    EmptyPSSMQuery,
    PSSMQueryTooLong {
        length: usize,
    },
}

impl Display for Error {
    fn fmt(&self, f: &mut Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{self:?}")
    }
}

impl std::error::Error for Error {}
