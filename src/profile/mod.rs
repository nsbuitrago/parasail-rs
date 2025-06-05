mod error;
mod profile_creator;

use crate::{InstructionSet, Matrix};
use libparasail_sys::{parasail_matrix_t, parasail_profile_t};
use profile_creator::profile_creator_lookup;
use std::ffi::{c_int, c_void, CString};
use std::ops::Deref;

pub use error::*; // flatten

/// Profile builder
#[derive(Debug, Clone)]
pub struct ProfileBuilder<'a> {
    query: &'a [u8],
    matrix: &'a Matrix,
    use_stats: bool,
    solution_width: Option<u8>,
    instruction_set: Option<InstructionSet>,
}

impl<'a> ProfileBuilder<'a> {
    /// Prepare profile for an alignment function that returns statistics.
    pub fn use_stats(&mut self) -> &mut Self {
        self.use_stats = true;
        self
    }

    /// Set solution width (i.e., 8. 16, 32, or 64 bits). If not set here, will
    /// allocate profiles for 8 and 16 bit solutions.
    pub fn solution_width(&mut self, width: u8) -> &mut Self {
        self.solution_width = Some(width);
        self
    }

    /// Explicitly define SIMD instruction set. If not set here, the best
    /// instruction set will be chosen automatically.
    pub fn instruction_set(&mut self, instruction_set: InstructionSet) -> &mut Self {
        self.instruction_set = Some(instruction_set);
        self
    }

    /// Build a new query profile.
    ///
    /// This method can fail if there is an invalid solution width,
    /// the query contains interior null bytes, profile creation function lookup fails or returns
    /// null pointer.
    pub fn build(&self) -> Result<Profile> {
        // The parasail_lookup_pcreator function does not work as expected for some
        // reason. We just use our own dispatcher here.
        //
        // let create_profile_fn = unsafe { parasail_lookup_pcreator(fn_name.as_ptr()) };
        // if create_profile_fn.is_none() {
        //     return Err(Error::FunctionLookupFailed {func_name: fn_name.to_string_lossy().into() });
        // }
        let profile_fn =
            profile_creator_lookup(self.use_stats, self.instruction_set, self.solution_width)?;

        let query_len = self.query.len() as c_int;
        let query_cstring = CString::new(self.query)?;
        let matrix = match self.matrix {
            Matrix::Builtin(matrix) => *matrix,
            Matrix::Custom(matrix) => *matrix as *const parasail_matrix_t,
        };

        let profile = unsafe { profile_fn(query_cstring.as_ptr(), query_len, matrix) };
        if profile.is_null() {
            Err(Error::NullProfile {
                msg: "Profile creation returned null pointer".to_string(),
            })
        } else {
            Ok(Profile {
                inner: profile,
                use_stats: self.use_stats,
            })
        }
    }
}

/// Query profile for alignment with striped or scan vectorization.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct Profile {
    pub(crate) inner: *mut parasail_profile_t,
    pub(crate) use_stats: bool,
}

impl Profile {
    /// Create a new profile builder.
    pub fn new<'a>(query: &'a [u8], matrix: &'a Matrix) -> ProfileBuilder<'a> {
        ProfileBuilder {
            query,
            matrix,
            use_stats: false,
            solution_width: None,
            instruction_set: None,
        }
    }
}

#[doc(hidden)]
impl Deref for Profile {
    type Target = *mut parasail_profile_t;
    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

#[doc(hidden)]
impl Drop for Profile {
    fn drop(&mut self) {
        unsafe {
            if let Some(free_fn) = (*self.inner).free {
                free_fn(self.inner as *mut c_void);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(target_arch = "aarch64")]
    #[test]
    fn function_lookup() {
        let profile_fn = profile_creator_lookup(false, Some(InstructionSet::AltiVec), Some(16));
        assert!(profile_fn.is_ok());
    }
}
