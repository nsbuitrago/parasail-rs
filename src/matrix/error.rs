use derive_more::From;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, From)]
pub enum Error {
    #[from]
    // An interior null byte is found in a string when converting to CString.
    NulError(std::ffi::NulError),
    // Built in matrix <name> not found.
    MatrixNotFound {
        matrix_name: String,
    },
    // Failed to lookup and create a builtin matrix.
    MatrixLookupFailed {
        matrix_name: String,
    },
    // Custom matrix creation failed.
    MatrixCreationFailed {
        alphabet: String,
        match_score: i32,
        mismatch_score: i32,
    },
    // PSSM creation failed
    PSSMCreationFailed {
        alphabet: String,
        values: Vec<i32>,
        n_rows: i32,
    },
    FromFileFailed {
        filepath: String,
    },
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{self:?}")
    }
}

impl std::error::Error for Error {}
