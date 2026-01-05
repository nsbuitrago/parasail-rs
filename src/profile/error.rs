use crate::{InstructionSet, SolutionWidth};
use derive_more::From;
use std::ffi::NulError;
use std::fmt::{Display, Formatter};

#[derive(Debug, From)]
pub enum Error {
    QueryIsEmpty,
    ProfileFnLookupFailed {
        use_stats: bool,
        instruction_set: InstructionSet,
        solution_width: SolutionWidth,
    },
    NullProfile,
    #[from]
    InteriorNulByte(NulError),
}

impl Display for Error {
    fn fmt(&self, f: &mut Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{self:?}")
    }
}

impl std::error::Error for Error {}
