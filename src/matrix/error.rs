use std::ffi::NulError;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum Error {
    #[error("Error creating matrix: {0}")]
    CreateErr(#[from] NulError),
    #[error("Error creating matrix from lookup: {0} matrix not found.")]
    LookupErr(String),
    #[error("Error creating matrix from file: {0}")]
    FromFileErr(String),
    #[error("File not found: {0}")]
    FileNotFound(String),
    #[error("Error creating PSSM matrix. Invalid alphabet, values, or rows.")]
    CreatePssmErr,
    #[error("Error creating matrix. Matrix has not been created yet. Consider using the `create` method.")]
    NullMatrix,
    #[error("Error converting matrix to PSSM. Matrix is not square.")]
    NotSquare,
    #[error("Error setting value on substitution matrix: Matrix must be a user matrix and not builtin. Consider using the `create` method")]
    NotBuiltIn,
    #[error("Error setting value on substitution matrix: Index ({0},{1}) out of range ((0,0),({2},{2}))")]
    InvalidIndex(i32, i32, i32),
}
