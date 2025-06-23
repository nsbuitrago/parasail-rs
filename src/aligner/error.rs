use std::ffi::NulError;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum Error {
    #[error("Alignment initialization error: {0}")]
    AlignInitErr(#[from] NulError),
    #[error("No bandwidth set for banded alignment.")]
    NoBandwidth,
    #[error("Alignment error: {0}")]
    AlignmentErr(#[from] crate::aligner::alignment::Error),
}
