pub mod error;

use crate::{InstructionSet, Matrix, Profile};
pub use error::*; // flatten
use libparasail_sys::{
    parasail_lookup_function, parasail_lookup_pfunction, parasail_matrix_t, parasail_profile_t,
    parasail_result_free, parasail_result_get_end_query, parasail_result_get_end_ref,
    parasail_result_get_length, parasail_result_get_matches, parasail_result_get_score,
    parasail_result_get_similar, parasail_result_is_banded, parasail_result_is_blocked,
    parasail_result_is_diag, parasail_result_is_nw, parasail_result_is_rowcol,
    parasail_result_is_saturated, parasail_result_is_scan, parasail_result_is_sg,
    parasail_result_is_stats, parasail_result_is_stats_rowcol, parasail_result_is_stats_table,
    parasail_result_is_striped, parasail_result_is_sw, parasail_result_is_table,
    parasail_result_is_trace, parasail_result_t,
};
use log::warn;
use std::ffi::{c_char, CString};
use std::sync::Arc;

/// A gap enum to specify where to allow or not penalize gaps.
pub enum Gap {
    Prefix,
    Suffix,
    All,
}

/// Aligner builder struct
#[non_exhaustive]
#[derive(Debug, Clone)]
pub struct AlignerBuilder {
    mode: &'static str,
    solution_width: Option<u8>,
    matrix: Arc<Matrix>,
    gap_open_penalty: u32,
    gap_extend_penalty: u32,
    profile: Arc<Option<Profile>>,
    allow_query_gaps: &'static str,
    allow_ref_gaps: &'static str,
    vec_strategy: &'static str,
    use_stats: &'static str,
    use_table: &'static str,
    use_trace: &'static str,
    instruction_set: &'static str,
}

impl AlignerBuilder {
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

    /// Set the scoring matrix used for alignment.
    pub fn matrix(&mut self, matrix: Matrix) -> &mut Self {
        self.matrix = Arc::new(matrix);
        self
    }

    /// Set a custom gap opening penalty for alignment. By default, gap openings
    /// are not penalized.
    pub fn gap_open_penalty(&mut self, penalty: u32) -> &mut Self {
        self.gap_open_penalty = penalty;
        self
    }

    /// Set a custom gap extension penalty for alignment. By default, Gap extensions
    /// are not penalized.
    pub fn gap_extend_penalty(&mut self, penalty: u32) -> &mut Self {
        self.gap_extend_penalty = penalty;
        self
    }

    /// Set which gaps are allowed (not penalized) on the query sequence for
    /// semi-global alignment. By default, gaps are penalized at the beginning
    /// and end (prefix and suffix) of the query sequence if penalties have been
    /// set.
    ///
    /// Example:
    /// ```rust, no_run
    /// use parasail_rs::{Aligner, Gap};
    /// // allow gaps on the prefix (i.e., beginning of the query)
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
    /// for the reference sequence.
    pub fn allow_ref_gaps(&mut self, gap: Gap) -> &mut Self {
        match gap {
            Gap::Prefix => self.allow_ref_gaps = "_db",
            Gap::Suffix => self.allow_ref_gaps = "_de",
            Gap::All => self.allow_ref_gaps = "_dx",
        }
        self
    }

    /// Return alignment statistics. Note that `return_stats` and `enable_traceback`
    /// are mutually exclusive. Calling this method will disable traceback
    /// if it has been enabled for this aligner previously.
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

    /// Enable traceback capable alignment. Note that `enable_traceback` and
    /// `return_stats` are mutually exclusive. Calling this method will disable
    /// returning statistics if it has been enabled for the aligner previously.
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

    /// Return the DP table. Note that if `return_last_rowcol` has been set previously,
    /// only the last row and column from the table will be return. Use `return_table`
    /// (without calling `return_last_rowcol`) to return the entire table.
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
        // prefer returning last row/col from DP table if set
        if self.use_table == "_table" {
            warn!(
                "Returning the entire table was previously enabled.
                To return the entire table, use `return_table` exclusively.
                Returning last row and col only."
            );
        }
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
    pub fn profile(&mut self, profile: Profile) -> &mut Self {
        // update use stats based on the profile
        if profile.use_stats {
            self.use_stats = "_stats";
        }

        self.profile = Arc::new(Some(profile));
        self
    }

    /// Create a new aligner.
    pub fn build(&self) -> Result<Aligner> {
        // get the function name for lookup
        let alignment_fn_name = self.get_parasail_fn_name()?;

        let parasail_fn = if self.profile.is_some() {
            // make sure the vectorization strategy is compatible with profile seting
            if self.vec_strategy != "_striped" && self.vec_strategy != "_scan" {
                println!("not recognizing vec strategy?, {:?}", self.vec_strategy);
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

        let aligner = Aligner {
            parasail_fn,
            matrix: Arc::clone(&self.matrix),
            profile: Arc::clone(&self.profile),
            gap_open_penalty: self.gap_open_penalty,
            gap_extend_penalty: self.gap_extend_penalty,
        };

        Ok(aligner)
    }

    /// Get the parasail function name from the set fields.
    fn get_parasail_fn_name(&self) -> Result<CString> {
        let (instruction_set, solution_width) = self.get_vectorization_config()?;

        let profile = if self.profile.is_some() {
            "_profile"
        } else {
            ""
        };

        let fn_name = CString::new(format!(
            "parasail_{}{}{}{}{}{}{}{}{}{}",
            self.mode,
            self.allow_query_gaps,
            self.allow_ref_gaps,
            self.use_stats,
            self.use_table,
            self.use_trace,
            self.vec_strategy,
            profile,
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

impl Default for AlignerBuilder {
    /// Default is an unvectorized global aligner with a DNA identity scoring
    /// matrix and no gap penalties.
    fn default() -> Self {
        AlignerBuilder {
            mode: "nw",
            solution_width: None,
            matrix: Arc::new(Matrix::default()),
            gap_open_penalty: 0,
            gap_extend_penalty: 0,
            profile: None.into(),
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
#[derive(Debug, Clone)]
pub struct Aligner {
    parasail_fn: AlignerFn,
    pub matrix: Arc<Matrix>,
    pub profile: Arc<Option<Profile>>,
    pub gap_open_penalty: u32,
    pub gap_extend_penalty: u32,
}

impl Aligner {
    /// Create a new aligner builder. Default is a global aligner with
    pub fn new() -> AlignerBuilder {
        AlignerBuilder::default()
    }

    /// Align query and reference sequences.
    pub fn align(&self, query: &[u8], reference: &[u8]) -> Result<Alignment> {
        let query_len = query.len();
        let reference_len = reference.len();
        let query = CString::new(query)?;
        let reference = CString::new(reference)?;

        let alignment = match self.parasail_fn {
            AlignerFn::Function(function) => {
                unsafe {
                    // we can safely unwrap since we checked the function is Some
                    // during the build call.
                    function.unwrap()(
                        query.as_ptr(),
                        query_len as i32,
                        reference.as_ptr(),
                        reference_len as i32,
                        self.gap_open_penalty as i32,
                        self.gap_extend_penalty as i32,
                        **self.matrix,
                    )
                }
            }
            AlignerFn::PFunction(_) => {
                return Err(Error::IncompatibleAlignerFn {
                    aligner_fn: self.parasail_fn,
                });
            }
        };

        Ok(Alignment { inner: alignment })
    }

    /// Align sequence with query profile
    pub fn align_with_profile(&self, reference: &[u8]) -> Result<Alignment> {
        let ref_seq_len = reference.len();
        let ref_seq = CString::new(reference)?;

        let alignment = match self.parasail_fn {
            AlignerFn::PFunction(function) => {
                let alignment = if let Some(ref profile) = *self.profile {
                    // let profile_ptr = profile.inner;
                    unsafe {
                        function.unwrap()(
                            profile.inner,
                            ref_seq.as_ptr(),
                            ref_seq_len as i32,
                            self.gap_open_penalty as i32,
                            self.gap_extend_penalty as i32,
                        )
                    }
                } else {
                    return Err(Error::NoProfileFound);
                };

                alignment
            }
            AlignerFn::Function(_) => {
                return Err(Error::IncompatibleAlignerFn {
                    aligner_fn: self.parasail_fn,
                });
            }
        };

        if alignment.is_null() {
            return Err(Error::NulResult {
                msg: "Alignment result is a null pointer".to_string(),
            });
        }

        Ok(Alignment { inner: alignment })
    }
}

#[doc(hidden)]
unsafe impl Send for Aligner {}
#[doc(hidden)]
unsafe impl Sync for Aligner {}

// Type of aligner function (either with or without profile)
#[derive(Debug, Clone, Copy)]
pub enum AlignerFn {
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

/// Alignment result returned from `align` or `align_with_profile` methods.
#[non_exhaustive]
#[derive(Debug)]
pub struct Alignment {
    inner: *const parasail_result_t,
}

impl Alignment {
    /// Get the alignment score
    pub fn get_score(&self) -> i32 {
        unsafe { parasail_result_get_score(self.inner) }
    }

    /// Get the zero-indexed position of the end of the query match
    pub fn get_end_query(&self) -> i32 {
        unsafe { parasail_result_get_end_query(self.inner) }
    }

    /// Get the zero-indexed position of the end of the reference match
    pub fn get_end_ref(&self) -> i32 {
        unsafe { parasail_result_get_end_ref(self.inner) }
    }

    /// Get the number of matches in the alignment.
    ///
    /// This method can return an error if the Aligner is set to not return statistics.
    /// The default behavior of the aligner is to not return statistics. To change
    /// this, call the `return_stats` method on the Aligner.
    pub fn get_matches(&self) -> Result<i32> {
        if self.is_stats() {
            unsafe { Ok(parasail_result_get_matches(self.inner)) }
        } else {
            Err(Error::NoStatsReturned)
        }
    }

    /// Get the number of similar residues (in the case of using a PSSM).
    ///
    /// This method can return an error if the Aligner is set to not return statistics.
    /// The default behavior of the aligner is to not return statistics. To change
    /// this, call the `return_stats` method on the Aligner.
    pub fn get_similar(&self) -> Result<i32> {
        if self.is_stats() {
            unsafe { Ok(parasail_result_get_similar(self.inner)) }
        } else {
            Err(Error::NoStatsReturned)
        }
    }

    /// Get alignment length
    ///
    /// This method can return an error if the Aligner is set to not return statistics.
    /// The default behavior of the aligner is to not return statistics. To change
    /// this, call the `return_stats` method on the Aligner.
    pub fn get_length(&self) -> Result<i32> {
        if self.is_stats() {
            unsafe { Ok(parasail_result_get_length(self.inner)) }
        } else {
            Err(Error::NoStatsReturned)
        }
    }

    /// Get score table
    pub fn get_score_table(&self) -> Result<i32> {
        let score_table = unsafe { parasail_result_get_score_table(self.inner) };
        return Ok(score_table);
    }

    pub fn get_matches_table(&self) -> Result<i32> {
        todo!()
    }

    pub fn get_similar_table(&self) -> Result<i32> {
        todo!()
    }

    pub fn get_length_table(&self) -> Result<i32> {
        todo!()
    }

    pub fn get_score_row(&self) -> Result<i32> {
        todo!()
    }

    pub fn get_matches_row(&self) -> Result<i32> {
        todo!()
    }

    pub fn get_similar_row(&self) -> Result<i32> {
        todo!()
    }

    pub fn get_length_row(&self) -> Result<i32> {
        todo!()
    }

    pub fn get_score_col(&self) -> Result<i32> {
        todo!()
    }

    pub fn get_matches_col(&self) -> Result<i32> {
        todo!()
    }

    pub fn get_similar_col(&self) -> Result<i32> {
        todo!()
    }

    pub fn get_length_col(&self) -> Result<i32> {
        todo!()
    }

    pub fn get_trace_table(&self) -> Result<i32> {
        todo!()
    }

    pub fn get_trace_ins_table(&self) -> Result<i32> {
        todo!()
    }

    pub fn get_trace_del_table(&self) -> Result<i32> {
        todo!()
    }

    /// This method can return an error if the Aligner is set to not return statistics.
    /// The default behavior of the aligner is to not return statistics. To change
    /// this, call the `return_stats` method on the Aligner.
    pub fn print_traceback() -> Result<()> {
        println!("get traceback string");
        Ok(())
    }

    pub fn get_traceback_strings() {
        todo!()
    }

    pub fn get_cigar() {
        todo!()
    }

    /// Check if the alignment mode is global.
    pub fn is_global(&self) -> bool {
        unsafe { parasail_result_is_nw(self.inner) != 0 }
    }

    /// Check if the alignment mode is semi-global.
    pub fn is_semi_global(&self) -> bool {
        unsafe { parasail_result_is_sg(self.inner) != 0 }
    }

    /// Check if the alignment mode is local.
    pub fn is_local(&self) -> bool {
        unsafe { parasail_result_is_sw(self.inner) != 0 }
    }

    /// Check if the solution width is saturated (i.e., using 8-bit solution width first and
    /// falling back to 16-bit if necessary).
    pub fn is_saturated(&self) -> bool {
        unsafe { parasail_result_is_saturated(self.inner) != 0 }
    }

    /// Check if banded alignment is used.
    pub fn is_banded(&self) -> bool {
        unsafe { parasail_result_is_banded(self.inner) != 0 }
    }

    /// Check if vector strategy is scan.
    pub fn is_scan(&self) -> bool {
        unsafe { parasail_result_is_scan(self.inner) != 0 }
    }

    /// Check if vector strategy is striped.
    pub fn is_striped(&self) -> bool {
        unsafe { parasail_result_is_striped(self.inner) != 0 }
    }

    /// Check if vector strategy is diagonal.
    pub fn is_diag(&self) -> bool {
        unsafe { parasail_result_is_diag(self.inner) != 0 }
    }

    pub fn is_blocked(&self) -> bool {
        unsafe { parasail_result_is_blocked(self.inner) != 0 }
    }

    /// Check if statistics are returned from alignment.
    pub fn is_stats(&self) -> bool {
        unsafe { parasail_result_is_stats(self.inner) != 0 }
    }

    /// Check if result is a stats table
    pub fn is_stats_table(&self) -> bool {
        unsafe { parasail_result_is_stats_table(self.inner) != 0 }
    }

    /// Check if result is a table
    pub fn is_table(&self) -> bool {
        unsafe { parasail_result_is_table(self.inner) != 0 }
    }

    /// Check if result is a last row and column of table
    pub fn is_rowcol(&self) -> bool {
        unsafe { parasail_result_is_rowcol(self.inner) != 0 }
    }

    /// Check if result is a row and column of table with additional statistics.
    pub fn is_stats_rowcol(&self) -> bool {
        unsafe { parasail_result_is_stats_rowcol(self.inner) != 0 }
    }

    /// Check if result is trace enabled.
    pub fn is_trace(&self) -> bool {
        unsafe { parasail_result_is_trace(self.inner) != 0 }
    }
}

#[doc(hidden)]
impl Drop for Alignment {
    fn drop(&mut self) {
        unsafe { parasail_result_free(self.inner as *mut parasail_result_t) }
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
