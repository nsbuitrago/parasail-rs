mod error;

use libparasail_sys::{
    parasail_matrix, parasail_profile, parasail_profile_create_16, parasail_profile_create_32,
    parasail_profile_create_64, parasail_profile_create_8, parasail_profile_create_altivec_128_16,
    parasail_profile_create_altivec_128_32, parasail_profile_create_altivec_128_64,
    parasail_profile_create_altivec_128_8, parasail_profile_create_altivec_128_sat,
    parasail_profile_create_avx_256_16, parasail_profile_create_avx_256_32,
    parasail_profile_create_avx_256_64, parasail_profile_create_avx_256_8,
    parasail_profile_create_avx_256_sat, parasail_profile_create_neon_128_16,
    parasail_profile_create_neon_128_32, parasail_profile_create_neon_128_64,
    parasail_profile_create_neon_128_8, parasail_profile_create_neon_128_sat,
    parasail_profile_create_sat, parasail_profile_create_sse_128_16,
    parasail_profile_create_sse_128_32, parasail_profile_create_sse_128_64,
    parasail_profile_create_sse_128_8, parasail_profile_create_sse_128_sat,
    parasail_profile_create_stats_16, parasail_profile_create_stats_32,
    parasail_profile_create_stats_64, parasail_profile_create_stats_8,
    parasail_profile_create_stats_altivec_128_16, parasail_profile_create_stats_altivec_128_32,
    parasail_profile_create_stats_altivec_128_64, parasail_profile_create_stats_altivec_128_8,
    parasail_profile_create_stats_altivec_128_sat, parasail_profile_create_stats_avx_256_16,
    parasail_profile_create_stats_avx_256_32, parasail_profile_create_stats_avx_256_64,
    parasail_profile_create_stats_avx_256_8, parasail_profile_create_stats_avx_256_sat,
    parasail_profile_create_stats_neon_128_16, parasail_profile_create_stats_neon_128_32,
    parasail_profile_create_stats_neon_128_64, parasail_profile_create_stats_neon_128_8,
    parasail_profile_create_stats_neon_128_sat, parasail_profile_create_stats_sat,
    parasail_profile_create_stats_sse_128_16, parasail_profile_create_stats_sse_128_32,
    parasail_profile_create_stats_sse_128_64, parasail_profile_create_stats_sse_128_8,
    parasail_profile_create_stats_sse_128_sat, parasail_profile_free, parasail_profile_t,
    parasail_ssw_init,
};
use std::ffi::{c_int, CString};
use std::ops::Deref;

use crate::{InstructionSet, Matrix, Result, SolutionWidth};

pub use error::*;

pub struct ProfileBuilder<'a> {
    query: &'a [u8],
    matrix: &'a Matrix,
    use_stats: bool,
    solution_width: SolutionWidth,
    instruction_set: InstructionSet,
}

impl<'a> ProfileBuilder<'a> {
    /// Create a new profile builder given a query and matrix. By default, the profile builder does
    /// not enable returning statistics and uses an 8-bit solution width with the option of trying
    /// 16-bit solution width in case of overflow. Additionally, the best SIMD instruction set is
    /// determined automatically. See additional methods for configuring the profile.
    pub fn new(query: &'a [u8], matrix: &'a Matrix) -> Self {
        Self {
            query,
            matrix,
            use_stats: false,
            solution_width: SolutionWidth::Sat,
            instruction_set: InstructionSet::Best,
        }
    }

    /// Prepare profile for an alignment that returns statistics
    pub fn use_stats(&mut self) -> &mut Self {
        self.use_stats = true;
        self
    }

    /// Set solution width (i.e., 8. 16, 32, or 64 bits). If not set here, will
    /// allocate profiles for 8 and 16 bit solutions.
    pub fn solution_width(&mut self, solution_width: SolutionWidth) -> &mut Self {
        self.solution_width = solution_width;
        self
    }

    /// Explicitly define SIMD instruction set. If not set here, the best
    /// instruction set will be chosen automatically.
    pub fn instruction_set(&mut self, instruction_set: InstructionSet) -> &mut Self {
        self.instruction_set = instruction_set;
        self
    }

    /// Build a new query profile.
    pub fn build(&self) -> Result<Profile> {
        let create_profile = self.profile_creator_lookup();

        let query_len = self.query.len() as c_int;
        let query_cstring = CString::new(self.query).map_err(Error::NulError)?;

        let profile = unsafe {
            create_profile(
                query_cstring.as_ptr(),
                self.query.len() as c_int,
                self.matrix.inner,
            )
        };

        if profile.is_null() {
            return Err(Error::NullProfile.into());
        }

        Ok(Profile {
            inner: profile,
            use_stats: self.use_stats,
            query_len,
        })
    }

    fn profile_creator_lookup(
        &self,
    ) -> unsafe extern "C" fn(*const i8, i32, *const parasail_matrix) -> *mut parasail_profile {
        match (&self.use_stats, &self.instruction_set, &self.solution_width) {
            // without stats
            (false, InstructionSet::SSE2, SolutionWidth::Sat) => {
                parasail_profile_create_sse_128_sat
            }
            (false, InstructionSet::SSE2, SolutionWidth::Bit8) => parasail_profile_create_sse_128_8,
            (false, InstructionSet::SSE2, SolutionWidth::Bit16) => {
                parasail_profile_create_sse_128_16
            }
            (false, InstructionSet::SSE2, SolutionWidth::Bit32) => {
                parasail_profile_create_sse_128_32
            }
            (false, InstructionSet::SSE2, SolutionWidth::Bit64) => {
                parasail_profile_create_sse_128_64
            }
            (false, InstructionSet::SSE41, SolutionWidth::Sat) => {
                parasail_profile_create_sse_128_sat
            }
            (false, InstructionSet::SSE41, SolutionWidth::Bit8) => {
                parasail_profile_create_sse_128_8
            }
            (false, InstructionSet::SSE41, SolutionWidth::Bit16) => {
                parasail_profile_create_sse_128_16
            }
            (false, InstructionSet::SSE41, SolutionWidth::Bit32) => {
                parasail_profile_create_sse_128_32
            }
            (false, InstructionSet::SSE41, SolutionWidth::Bit64) => {
                parasail_profile_create_sse_128_64
            }
            (false, InstructionSet::AVX2, SolutionWidth::Sat) => {
                parasail_profile_create_avx_256_sat
            }
            (false, InstructionSet::AVX2, SolutionWidth::Bit8) => parasail_profile_create_avx_256_8,
            (false, InstructionSet::AVX2, SolutionWidth::Bit16) => {
                parasail_profile_create_avx_256_16
            }
            (false, InstructionSet::AVX2, SolutionWidth::Bit32) => {
                parasail_profile_create_avx_256_32
            }
            (false, InstructionSet::AVX2, SolutionWidth::Bit64) => {
                parasail_profile_create_avx_256_64
            }
            (false, InstructionSet::Neon, SolutionWidth::Sat) => {
                parasail_profile_create_neon_128_sat
            }
            (false, InstructionSet::Neon, SolutionWidth::Bit8) => {
                parasail_profile_create_neon_128_8
            }
            (false, InstructionSet::Neon, SolutionWidth::Bit16) => {
                parasail_profile_create_neon_128_16
            }
            (false, InstructionSet::Neon, SolutionWidth::Bit32) => {
                parasail_profile_create_neon_128_32
            }
            (false, InstructionSet::Neon, SolutionWidth::Bit64) => {
                parasail_profile_create_neon_128_64
            }
            (false, InstructionSet::AltiVec, SolutionWidth::Sat) => {
                parasail_profile_create_altivec_128_sat
            }
            (false, InstructionSet::AltiVec, SolutionWidth::Bit8) => {
                parasail_profile_create_altivec_128_8
            }
            (false, InstructionSet::AltiVec, SolutionWidth::Bit16) => {
                parasail_profile_create_altivec_128_16
            }
            (false, InstructionSet::AltiVec, SolutionWidth::Bit32) => {
                parasail_profile_create_altivec_128_32
            }
            (false, InstructionSet::AltiVec, SolutionWidth::Bit64) => {
                parasail_profile_create_altivec_128_64
            }
            (false, InstructionSet::Best, SolutionWidth::Sat) => parasail_profile_create_sat,
            (false, InstructionSet::Best, SolutionWidth::Bit8) => parasail_profile_create_8,
            (false, InstructionSet::Best, SolutionWidth::Bit16) => parasail_profile_create_16,
            (false, InstructionSet::Best, SolutionWidth::Bit32) => parasail_profile_create_32,
            (false, InstructionSet::Best, SolutionWidth::Bit64) => parasail_profile_create_64,
            // with stats
            (true, InstructionSet::SSE2, SolutionWidth::Sat) => {
                parasail_profile_create_stats_sse_128_sat
            }
            (true, InstructionSet::SSE2, SolutionWidth::Bit8) => {
                parasail_profile_create_stats_sse_128_8
            }
            (true, InstructionSet::SSE2, SolutionWidth::Bit16) => {
                parasail_profile_create_stats_sse_128_16
            }
            (true, InstructionSet::SSE2, SolutionWidth::Bit32) => {
                parasail_profile_create_stats_sse_128_32
            }
            (true, InstructionSet::SSE2, SolutionWidth::Bit64) => {
                parasail_profile_create_stats_sse_128_64
            }
            (true, InstructionSet::SSE41, SolutionWidth::Sat) => {
                parasail_profile_create_stats_sse_128_sat
            }
            (true, InstructionSet::SSE41, SolutionWidth::Bit8) => {
                parasail_profile_create_stats_sse_128_8
            }
            (true, InstructionSet::SSE41, SolutionWidth::Bit16) => {
                parasail_profile_create_stats_sse_128_16
            }
            (true, InstructionSet::SSE41, SolutionWidth::Bit32) => {
                parasail_profile_create_stats_sse_128_32
            }
            (true, InstructionSet::SSE41, SolutionWidth::Bit64) => {
                parasail_profile_create_stats_sse_128_64
            }
            (true, InstructionSet::AVX2, SolutionWidth::Sat) => {
                parasail_profile_create_stats_avx_256_sat
            }
            (true, InstructionSet::AVX2, SolutionWidth::Bit8) => {
                parasail_profile_create_stats_avx_256_8
            }
            (true, InstructionSet::AVX2, SolutionWidth::Bit16) => {
                parasail_profile_create_stats_avx_256_16
            }
            (true, InstructionSet::AVX2, SolutionWidth::Bit32) => {
                parasail_profile_create_stats_avx_256_32
            }
            (true, InstructionSet::AVX2, SolutionWidth::Bit64) => {
                parasail_profile_create_stats_avx_256_64
            }
            (true, InstructionSet::Neon, SolutionWidth::Sat) => {
                parasail_profile_create_stats_neon_128_sat
            }
            (true, InstructionSet::Neon, SolutionWidth::Bit8) => {
                parasail_profile_create_stats_neon_128_8
            }
            (true, InstructionSet::Neon, SolutionWidth::Bit16) => {
                parasail_profile_create_stats_neon_128_16
            }
            (true, InstructionSet::Neon, SolutionWidth::Bit32) => {
                parasail_profile_create_stats_neon_128_32
            }
            (true, InstructionSet::Neon, SolutionWidth::Bit64) => {
                parasail_profile_create_stats_neon_128_64
            }
            (true, InstructionSet::AltiVec, SolutionWidth::Sat) => {
                parasail_profile_create_stats_altivec_128_sat
            }
            (true, InstructionSet::AltiVec, SolutionWidth::Bit8) => {
                parasail_profile_create_stats_altivec_128_8
            }
            (true, InstructionSet::AltiVec, SolutionWidth::Bit16) => {
                parasail_profile_create_stats_altivec_128_16
            }
            (true, InstructionSet::AltiVec, SolutionWidth::Bit32) => {
                parasail_profile_create_stats_altivec_128_32
            }
            (true, InstructionSet::AltiVec, SolutionWidth::Bit64) => {
                parasail_profile_create_stats_altivec_128_64
            }
            (true, InstructionSet::Best, SolutionWidth::Sat) => parasail_profile_create_stats_sat,
            (true, InstructionSet::Best, SolutionWidth::Bit8) => parasail_profile_create_stats_8,
            (true, InstructionSet::Best, SolutionWidth::Bit16) => parasail_profile_create_stats_16,
            (true, InstructionSet::Best, SolutionWidth::Bit32) => parasail_profile_create_stats_32,
            (true, InstructionSet::Best, SolutionWidth::Bit64) => parasail_profile_create_stats_64,
        }
    }
}

/// Query profile for sequence alignment
pub struct Profile {
    pub(crate) inner: *mut parasail_profile_t,
    pub(crate) use_stats: bool,
    pub(crate) query_len: i32,
}

impl Profile {
    pub fn builder<'a>(query: &'a [u8], matrix: &'a Matrix) -> ProfileBuilder<'a> {
        ProfileBuilder::new(query, matrix)
    }
    /// Create a new profile from a query sequence, to use with or without stats, and a scoring matrix.
    /// The with_stats should be set to true if you will use an alignment function that returns
    /// statistics. If true, the Profile will use the appropriate parasail functions to allocate
    /// additional data structures required for statistics.
    pub fn new(query_bytes: &[u8], with_stats: bool, matrix: &Matrix) -> Result<Self> {
        if query_bytes.is_empty() {
            return Err(Error::QueryIsEmpty.into());
        }

        let query_len = query_bytes.len() as i32;
        let query = CString::new(query_bytes).map_err(Error::NulError)?;

        unsafe {
            match with_stats {
                true => {
                    let profile =
                        parasail_profile_create_stats_sat(query.as_ptr(), query_len, **matrix);
                    if profile.is_null() {
                        return Err(Error::NullProfile.into());
                    }

                    Ok(Profile {
                        inner: profile,
                        use_stats: true,
                        query_len,
                    })
                }
                false => {
                    let profile = parasail_profile_create_sat(query.as_ptr(), query_len, **matrix);
                    if profile.is_null() {
                        return Err(Error::NullProfile.into());
                    }

                    Ok(Profile {
                        inner: profile,
                        use_stats: false,
                        query_len,
                    })
                }
            }
        }
    }

    pub fn new_ssw(query_bytes: &[u8], matrix: &Matrix, score_size: i8) -> Result<Self> {
        let query_len = query_bytes.len() as i32;
        if query_len == 0 {
            panic!("Query sequence has length 0.");
        }
        let query = CString::new(query_bytes).map_err(Error::NulError)?;

        let profile = unsafe {
            let profile = parasail_ssw_init(query.as_ptr(), query_len, **matrix, score_size);

            if profile.is_null() {
                return Err(Error::NullProfile.into());
            }
            profile
        };

        Ok(Profile {
            inner: profile,
            use_stats: true,
            query_len,
        })
    }
}

/// Default profile is a null pointer
// This is for cases where the profile is not used by the Aligner, but we need
// some default anyway. We could probably also not have this default and
// just wrap the Profile in an Option
impl Default for Profile {
    fn default() -> Self {
        Profile {
            inner: std::ptr::null_mut(),
            use_stats: false,
            query_len: 0,
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
        if !self.inner.is_null() {
            unsafe { parasail_profile_free(self.inner) }
        }
    }
}

#[doc(hidden)]
unsafe impl Send for Profile {}
#[doc(hidden)]
unsafe impl Sync for Profile {}
