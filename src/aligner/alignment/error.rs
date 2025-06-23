use std::ffi::{IntoStringError, NulError};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum Error {
    #[error("Error getting result from {0}. Stats must be initialized. Consider using `use_stats` method on AlignerBuilder.")]
    NoStats(String),
    #[error("Error getting result from {0}. Table must be enabled. Consider using `use_table` method on AlignerBuilder.")]
    NoTable(String),
    #[error("Error getting result from {0}. Table and stats must be enabled. Consider chaining `use_stats` and use_table methods on AlignerBuilder.")]
    NoStatsTable(String),
    #[error("Error getting result from {0}. Last row and col must be enabled. Consider using `use_rowcol` method on AlignerBuilder.")]
    NoRowCol(String),
    #[error("Error getting result from {0}. Traceback must be enabled. Consider using `use_trace` method on AlignerBuilder.")]
    NoTrace(String),
    #[error("Error converting CIGAR string to Rust string: {0}")]
    CigarToStringErr(#[from] IntoStringError),
    #[error("Error creating new CString: {0}")]
    NewCStringErr(#[from] NulError),
    #[error("No bandwidth set for banded alignment.")]
    NoBandwidth,
}
