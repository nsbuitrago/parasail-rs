pub mod error;

use crate::{InstructionSet, Matrix, Profile};
pub use error::*; // flatten
use libparasail_sys::{
    parasail_lookup_function, parasail_lookup_pfunction, parasail_matrix_t, parasail_profile_t,
    parasail_result_t,
};
use log::warn;
use std::ffi::{c_char, CString};

/// A gap enum to specify where to allow or not penalize gaps
pub enum Gap {
    Prefix,
    Suffix,
    All,
}

#[non_exhaustive]
#[derive(Debug, Clone)]
pub struct AlignerBuilder<'a> {
    mode: &'static str,
    solution_width: Option<u8>,
    matrix: Matrix,
    gap_open_penalty: i32,
    gap_extend_penalty: i32,
    profile: Option<&'a Profile>,
    allow_query_gaps: &'static str,
    allow_ref_gaps: &'static str,
    vec_strategy: &'static str,
    use_stats: &'static str,
    use_table: &'static str,
    use_trace: &'static str,
    instruction_set: &'static str,
}

impl<'a> AlignerBuilder<'a> {
    /// Use global (Needleman-Wunsch) alignment mode.
    pub fn global(&mut self) -> &mut Self {
        self.mode = "nw";
        self
    }

    /// Use semi-global alignment mode.
    pub fn semi_global(&mut self) -> &mut Self {
        self.mode = "sg";
        self
    }

    /// Use local (Smith-Watermann) alignment mode.
    pub fn local(&mut self) -> &mut Self {
        self.mode = "sw";
        self
    }

    /// Set a custom scoring matrix for alignment.
    pub fn matrix(&mut self, matrix: Matrix) -> &mut Self {
        self.matrix = matrix;
        self
    }

    /// Set a custom gap opening penalty for alignment. Gap openings are not
    /// penalized by default.
    pub fn gap_open_penalty(&mut self, penalty: i32) -> &mut Self {
        self.gap_open_penalty = penalty;
        self
    }

    /// Set a custom gap extension penalty for alignment. Gap extensions are not
    /// penalized by default.
    pub fn gap_extend_penalty(&mut self, penalty: i32) -> &mut Self {
        self.gap_extend_penalty = penalty;
        self
    }

    /// Set which gaps are allowed and not penalized on the query sequence for
    /// semi-global alignment. Gaps are penalized at the beginning and end
    /// (prefix and suffix) of the query sequence if penalties have been set by
    /// default. For example, to not penalize gaps at the beginning of the query
    /// sequence:
    /// Example:
    /// ```rust, no_run
    /// use parasail_rs::{Aligner, Gap};
    /// // to allow gaps on the prefix (i.e., beginning of the query)
    /// let aligner = Aligner::new().allow_query_gaps(Gap::Prefix).build()?;
    /// # Ok::<(), parasail_rs::Error>(())
    ///
    /// ```
    pub fn allow_query_gaps(&mut self, gap: Gap) -> &mut Self {
        match gap {
            Gap::Prefix => self.allow_query_gaps = "_qb",
            Gap::Suffix => self.allow_query_gaps = "_qe",
            Gap::All => self.allow_query_gaps = "_qx",
        }
        self
    }

    /// Set which gaps are allowed and not penalized on the reference sequence
    /// for semi-global alignment. Behaves the same as `allow_query_gaps`, but
    /// for the reference.
    pub fn allow_ref_gaps(&mut self, gap: Gap) -> &mut Self {
        match gap {
            Gap::Prefix => self.allow_ref_gaps = "_db",
            Gap::Suffix => self.allow_ref_gaps = "_de",
            Gap::All => self.allow_ref_gaps = "_dx",
        }
        self
    }

    /// Return alignment statistics.
    pub fn return_stats(&mut self) -> &mut Self {
        self.use_stats = "_stats";

        // disable traceback if already enabled
        if !self.use_trace.is_empty() {
            warn!(
                "Traceback was previously enabled, but mutually exclusive with
	        stats. Disabling traceback in favor of stats."
            );
            self.use_trace = "";
        }

        self
    }

    /// Enable traceback capable alignment
    pub fn enable_traceback(&mut self) -> &mut Self {
        self.use_trace = "_trace";

        // disable statistics if already enabled
        if !self.use_stats.is_empty() {
            warn!(
                "Alignment stats were previously enabled, but mutually exclusive
                with traceback. Disabling stats in favor of traceback."
            );
            self.use_stats = "";
        }
        self
    }

    /// Return the DP table
    pub fn return_table(&mut self) -> &mut Self {
        // prefer returning last row/col from DP table if set
        if self.use_table == "_rowcol" {
            warn!(
                "Returning last row/col from table was previously enabled.
                To return the entire table, use `return_table` exclusively.
                Returning last row and col only."
            );
        } else {
            self.use_table = "_table";
        }
        self
    }

    /// Return the last row/column of the DP table.
    pub fn return_last_rowcol(&mut self) -> &mut Self {
        self.use_table = "_rowcol";
        self
    }

    /// Use scan vectorization strategy.
    pub fn scan(&mut self) -> &mut Self {
        self.vec_strategy = "_scan";
        self
    }

    /// Use striped vectorization strategy.
    pub fn striped(&mut self) -> &mut Self {
        self.vec_strategy = "_striped";
        self
    }

    /// Use diag vectorization strategy.
    pub fn diag(&mut self) -> &mut Self {
        self.vec_strategy = "_diag";
        self
    }

    /// Set solution width (i.e., 8, 16, 32, 64).
    pub fn solution_width(&mut self, width: u8) -> &mut Self {
        self.solution_width = Some(width);
        self
    }

    pub fn instruction_set(&mut self, instruction_set: InstructionSet) -> &mut Self {
        match instruction_set {
            InstructionSet::SSE2 => self.instruction_set = "_sse2_128",
            InstructionSet::SSE41 => self.instruction_set = "_sse41_128",
            InstructionSet::AVX2 => self.instruction_set = "_avx2_256",
            InstructionSet::AltiVec => self.instruction_set = "_altivec_128",
            InstructionSet::Neon => self.instruction_set = "_neon_128",
        }
        self
    }

    /// Set query profile for alignment.
    pub fn profile(&mut self, profile: &'a Profile) -> &mut Self {
        self.profile = Some(profile);

        // update use stats based on the profile
        if profile.use_stats {
            self.use_stats = "_stats";
        }

        self
    }

    pub fn build(&self) -> Result<()> {
        // get the function name for lookup
        let alignment_fn_name = self.get_parasail_fn_name()?;

        let parasail_fn = if self.profile.is_some() {
            // make sure the vectorization strategy is compatible with profile seting
            if self.vec_strategy != "_striped" || self.vec_strategy != "_scan" {
                return Err(Error::IncompatibleProfileConfig {
                    vec_strategy: self.vec_strategy.to_string(),
                });
            }
            unsafe { AlignerFn::PFunction(parasail_lookup_pfunction(alignment_fn_name.as_ptr())) }
        } else {
            unsafe { AlignerFn::Function(parasail_lookup_function(alignment_fn_name.as_ptr())) }
        };

        if parasail_fn.is_none() {
            return Err(Error::FunctionLookupFailed {
                fn_name: alignment_fn_name,
            });
        }

        Ok(())
    }

    /// Get the parasail function name from the set fields.
    fn get_parasail_fn_name(&self) -> Result<CString> {
        let (instruction_set, solution_width) = self.get_vectorization_config()?;
        let fn_name = CString::new(format!(
            "parasail_{}{}{}{}{}{}{}{}{}",
            self.mode,
            self.allow_query_gaps,
            self.allow_ref_gaps,
            self.use_stats,
            self.use_table,
            self.use_trace,
            self.vec_strategy,
            instruction_set,
            solution_width
        ))?;

        Ok(fn_name)
    }

    /// Get the vectorization config string
    fn get_vectorization_config(&self) -> Result<(String, String)> {
        if self.vec_strategy.is_empty() {
            Ok((String::default(), String::default()))
        } else {
            let solution_width = match self.solution_width {
                Some(8 | 16 | 32 | 64) => {
                    let width_str = self.solution_width.unwrap().to_string();
                    format!("_{width_str}")
                }
                None => "_sat".to_string(),
                _ => {
                    return Err(Error::InvalidSolutionWidth {
                        solution_width: self.solution_width.unwrap(),
                    })
                }
            };
            Ok((self.instruction_set.to_string(), solution_width))
        }
    }
}

impl<'a> Default for AlignerBuilder<'a> {
    /// Default is an unvectorized global aligner with a DNA identity scoring
    /// matrix and no gap penalties.
    fn default() -> Self {
        AlignerBuilder {
            mode: "nw",
            solution_width: None,
            matrix: Matrix::default(),
            gap_open_penalty: 0,
            gap_extend_penalty: 0,
            profile: None,
            allow_query_gaps: "",
            allow_ref_gaps: "",
            vec_strategy: "",
            use_stats: "",
            use_table: "",
            use_trace: "",
            instruction_set: "",
        }
    }
}

/// Pairwise sequence aligner.
pub struct Aligner<'a> {
    parasail_fn: AlignerFn,
    pub matrix: &'a Matrix,
    pub profile: &'a Profile,
    pub gap_open_penalty: u32,
    pub gap_extend_penalty: u32,
}

impl<'a> Aligner<'a> {
    /// Create a new aligner builder. Default is a global aligner with
    pub fn new() -> AlignerBuilder<'a> {
        AlignerBuilder::default()
    }

    /// Align query and reference sequences.
    pub fn align(&self, query: &[u8], reference: &[u8]) -> Result<()> {
        todo!()
    }

    /// Align sequence with query profile
    pub fn align_with_profile(&self, reference: &[u8]) -> Result<()> {
        todo!()
    }
}

// Type of aligner function (either with or without profile)
enum AlignerFn {
    Function(
        Option<
            unsafe extern "C" fn(
                *const c_char,
                i32,
                *const c_char,
                i32,
                i32,
                i32,
                *const parasail_matrix_t,
            ) -> *mut parasail_result_t,
        >,
    ),
    PFunction(
        Option<
            unsafe extern "C" fn(
                *const parasail_profile_t,
                *const c_char,
                i32,
                i32,
                i32,
            ) -> *mut parasail_result_t,
        >,
    ),
}

/// Check if the aligner function is none.
impl AlignerFn {
    fn is_none(&self) -> bool {
        match self {
            AlignerFn::Function(f) => f.is_none(),
            AlignerFn::PFunction(f) => f.is_none(),
        }
    }
}

/// Ok alignment result returned from `align`  or `align_with_profile`.
#[non_exhaustive]
#[derive(Debug, Clone, Copy)]
pub struct Alignment {
    inner: f64,
}

impl Alignment {
    /// Get the alignment score
    pub fn score(&self) -> f64 {
        todo!();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn get_parasail_fn_name() -> Result<()> {
        // test basic configurations and then a few tricky ones

        // default NW aligner
        assert_eq!(
            Aligner::new().get_parasail_fn_name()?,
            CString::new("parasail_nw")?
        );

        // use local alignment
        assert_eq!(
            Aligner::new().local().get_parasail_fn_name()?,
            CString::new("parasail_sw")?
        );

        // use semi-global alignment
        assert_eq!(
            Aligner::new().semi_global().get_parasail_fn_name()?,
            CString::new("parasail_sg")?
        );

        // semi-global with allowed gaps
        assert_eq!(
            Aligner::new()
                .semi_global()
                .allow_query_gaps(Gap::Prefix)
                .get_parasail_fn_name()?,
            CString::new("parasail_sg_qb")?
        );

        assert_eq!(
            Aligner::new()
                .semi_global()
                .allow_ref_gaps(Gap::Suffix)
                .get_parasail_fn_name()?,
            CString::new("parasail_sg_de")?
        );

        // check that overriding table works when calling row_col
        assert_eq!(
            Aligner::new()
                .local()
                .return_table() // initially set to return full table
                .return_last_rowcol() // override and return the last row/col
                .get_parasail_fn_name()?,
            CString::new("parasail_sw_rowcol")?
        );

        // does this work if we do it the other way now,
        assert_eq!(
            Aligner::new()
                .local()
                .return_last_rowcol() // initially set to return last row/col
                .return_table() // call to return whole table but prefers row/col
                .get_parasail_fn_name()?,
            CString::new("parasail_sw_rowcol")?
        );

        // set a solution width along with some instructions
        // this should still work because I check for the vectorization strategy
        // first, and if it is not set, then we ignore the instruction set
        // and solution width settings.
        assert_eq!(
            Aligner::new()
                .global()
                .solution_width(8)
                .instruction_set(InstructionSet::Neon)
                .get_parasail_fn_name()?,
            CString::new("parasail_nw")?
        );

        // now we add the vec strategy and it should add the instruction set
        // and solution width information
        assert_eq!(
            Aligner::new()
                .global()
                .striped()
                .solution_width(8)
                .instruction_set(InstructionSet::Neon)
                .get_parasail_fn_name()?,
            CString::new("parasail_nw_striped_neon_128_8")?
        );

        Ok(())
    }
}
