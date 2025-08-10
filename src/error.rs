use derive_more::From;
pub type Result<T> = std::result::Result<T, Error>;

/// Parasail-rs errors.
#[derive(Debug, From)]
pub enum Error {
    #[from]
    Matrix(crate::matrix::Error),
    #[from]
    Profile(crate::profile::Error),
    #[from]
    Aligner(crate::aligner::Error),
    #[from]
    Alignment(crate::aligner::alignment::Error),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{self:?}")
    }
}

impl std::error::Error for Error {}
