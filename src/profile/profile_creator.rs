use crate::profile::{Error, Result};
use crate::InstructionSet;
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
    parasail_profile_create_stats_sse_128_sat,
};

pub fn profile_creator_lookup(
    use_stats: bool,
    instruction_set: Option<InstructionSet>,
    solution_width: Option<u8>,
) -> Result<unsafe extern "C" fn(*const i8, i32, *const parasail_matrix) -> *mut parasail_profile> {
    match (use_stats, instruction_set, solution_width) {
        (false, Some(InstructionSet::SSE2), None) => Ok(parasail_profile_create_sse_128_sat),
        (false, Some(InstructionSet::SSE2), Some(8)) => Ok(parasail_profile_create_sse_128_8),
        (false, Some(InstructionSet::SSE2), Some(16)) => Ok(parasail_profile_create_sse_128_16),
        (false, Some(InstructionSet::SSE2), Some(32)) => Ok(parasail_profile_create_sse_128_32),
        (false, Some(InstructionSet::SSE2), Some(64)) => Ok(parasail_profile_create_sse_128_64),
        (false, Some(InstructionSet::SSE41), None) => Ok(parasail_profile_create_sse_128_sat),
        (false, Some(InstructionSet::SSE41), Some(8)) => Ok(parasail_profile_create_sse_128_8),
        (false, Some(InstructionSet::SSE41), Some(16)) => Ok(parasail_profile_create_sse_128_16),
        (false, Some(InstructionSet::SSE41), Some(32)) => Ok(parasail_profile_create_sse_128_32),
        (false, Some(InstructionSet::SSE41), Some(64)) => Ok(parasail_profile_create_sse_128_64),
        (false, Some(InstructionSet::AVX2), None) => Ok(parasail_profile_create_avx_256_sat),
        (false, Some(InstructionSet::AVX2), Some(8)) => Ok(parasail_profile_create_avx_256_8),
        (false, Some(InstructionSet::AVX2), Some(16)) => Ok(parasail_profile_create_avx_256_16),
        (false, Some(InstructionSet::AVX2), Some(32)) => Ok(parasail_profile_create_avx_256_32),
        (false, Some(InstructionSet::AVX2), Some(64)) => Ok(parasail_profile_create_avx_256_64),
        (false, Some(InstructionSet::Neon), None) => Ok(parasail_profile_create_neon_128_sat),
        (false, Some(InstructionSet::Neon), Some(8)) => Ok(parasail_profile_create_neon_128_8),
        (false, Some(InstructionSet::Neon), Some(16)) => Ok(parasail_profile_create_neon_128_16),
        (false, Some(InstructionSet::Neon), Some(32)) => Ok(parasail_profile_create_neon_128_32),
        (false, Some(InstructionSet::Neon), Some(64)) => Ok(parasail_profile_create_neon_128_64),
        (false, Some(InstructionSet::AltiVec), None) => Ok(parasail_profile_create_altivec_128_sat),
        (false, Some(InstructionSet::AltiVec), Some(8)) => {
            Ok(parasail_profile_create_altivec_128_8)
        }

        (false, Some(InstructionSet::AltiVec), Some(16)) => {
            Ok(parasail_profile_create_altivec_128_16)
        }
        (false, Some(InstructionSet::AltiVec), Some(32)) => {
            Ok(parasail_profile_create_altivec_128_32)
        }
        (false, Some(InstructionSet::AltiVec), Some(64)) => {
            Ok(parasail_profile_create_altivec_128_64)
        }
        (false, None, None) => Ok(parasail_profile_create_sat),
        (false, None, Some(8)) => Ok(parasail_profile_create_8),
        (false, None, Some(16)) => Ok(parasail_profile_create_16),
        (false, None, Some(32)) => Ok(parasail_profile_create_32),
        (false, None, Some(64)) => Ok(parasail_profile_create_64),
        (true, Some(InstructionSet::SSE2), None) => Ok(parasail_profile_create_stats_sse_128_sat),
        (true, Some(InstructionSet::SSE2), Some(8)) => Ok(parasail_profile_create_stats_sse_128_8),
        (true, Some(InstructionSet::SSE2), Some(16)) => {
            Ok(parasail_profile_create_stats_sse_128_16)
        }
        (true, Some(InstructionSet::SSE2), Some(32)) => {
            Ok(parasail_profile_create_stats_sse_128_32)
        }
        (true, Some(InstructionSet::SSE2), Some(64)) => {
            Ok(parasail_profile_create_stats_sse_128_64)
        }
        (true, Some(InstructionSet::SSE41), None) => Ok(parasail_profile_create_stats_sse_128_sat),
        (true, Some(InstructionSet::SSE41), Some(8)) => Ok(parasail_profile_create_stats_sse_128_8),
        (true, Some(InstructionSet::SSE41), Some(16)) => {
            Ok(parasail_profile_create_stats_sse_128_16)
        }
        (true, Some(InstructionSet::SSE41), Some(32)) => {
            Ok(parasail_profile_create_stats_sse_128_32)
        }
        (true, Some(InstructionSet::SSE41), Some(64)) => {
            Ok(parasail_profile_create_stats_sse_128_64)
        }
        (true, Some(InstructionSet::AVX2), None) => Ok(parasail_profile_create_stats_avx_256_sat),
        (true, Some(InstructionSet::AVX2), Some(8)) => Ok(parasail_profile_create_stats_avx_256_8),
        (true, Some(InstructionSet::AVX2), Some(16)) => {
            Ok(parasail_profile_create_stats_avx_256_16)
        }
        (true, Some(InstructionSet::AVX2), Some(32)) => {
            Ok(parasail_profile_create_stats_avx_256_32)
        }
        (true, Some(InstructionSet::AVX2), Some(64)) => {
            Ok(parasail_profile_create_stats_avx_256_64)
        }
        (true, Some(InstructionSet::Neon), None) => Ok(parasail_profile_create_stats_neon_128_sat),
        (true, Some(InstructionSet::Neon), Some(8)) => Ok(parasail_profile_create_stats_neon_128_8),
        (true, Some(InstructionSet::Neon), Some(16)) => {
            Ok(parasail_profile_create_stats_neon_128_16)
        }
        (true, Some(InstructionSet::Neon), Some(32)) => {
            Ok(parasail_profile_create_stats_neon_128_32)
        }
        (true, Some(InstructionSet::Neon), Some(64)) => {
            Ok(parasail_profile_create_stats_neon_128_64)
        }
        (true, Some(InstructionSet::AltiVec), None) => {
            Ok(parasail_profile_create_stats_altivec_128_sat)
        }
        (true, Some(InstructionSet::AltiVec), Some(8)) => {
            Ok(parasail_profile_create_stats_altivec_128_8)
        }
        (true, Some(InstructionSet::AltiVec), Some(16)) => {
            Ok(parasail_profile_create_stats_altivec_128_16)
        }
        (true, Some(InstructionSet::AltiVec), Some(32)) => {
            Ok(parasail_profile_create_stats_altivec_128_32)
        }
        (true, Some(InstructionSet::AltiVec), Some(64)) => {
            Ok(parasail_profile_create_stats_altivec_128_64)
        }
        (true, None, None) => Ok(parasail_profile_create_stats_sat),
        (true, None, Some(8)) => Ok(parasail_profile_create_stats_8),
        (true, None, Some(16)) => Ok(parasail_profile_create_stats_16),
        (true, None, Some(32)) => Ok(parasail_profile_create_stats_32),
        (true, None, Some(64)) => Ok(parasail_profile_create_stats_64),
        _ => Err(Error::FunctionLookupFailed {
            use_stats,
            instruction_set,
            solution_width,
        }),
    }
}
