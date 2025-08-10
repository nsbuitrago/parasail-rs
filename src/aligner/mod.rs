use crate::{Matrix, Profile};
use libc::c_char;
use libparasail_sys::{
    parasail_lookup_function, parasail_lookup_pfunction, parasail_matrix_t, parasail_nw_banded,
    parasail_profile_t, parasail_result_t, parasail_ssw,
};
use log::warn;
use std::ffi::CString;
use std::sync::Arc;

pub(crate) mod alignment;
mod error;

use crate::Result;
pub use alignment::*;
pub use error::Error;

/// Parasail alignment function type.
#[derive(Clone)]
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

/// Aligner builder
pub struct AlignerBuilder {
    mode: String,
    solution_width: String,
    matrix: Arc<Matrix>,
    gap_open: i32,
    gap_extend: i32,
    profile: Arc<Profile>,
    allow_query_gaps: Vec<String>,
    allow_ref_gaps: Vec<String>,
    vec_strategy: String,
    use_stats: String,
    use_table: String,
    use_trace: String,
    bandwidth: Option<i32>,
}

/// Default aligner uses global alignment with an identity matrix for DNA
/// sequences and no gap penalties. No profile, trace, table, or stats options
/// are set. Vectorization strategy is set to striped by default.
impl Default for AlignerBuilder {
    fn default() -> Self {
        AlignerBuilder {
            mode: String::from("nw"),
            solution_width: String::from("sat"),
            matrix: Matrix::default().into(),
            gap_open: 0,
            gap_extend: 0,
            profile: Profile::default().into(),
            allow_query_gaps: Vec::default(),
            allow_ref_gaps: Vec::default(),
            vec_strategy: String::from("_striped"),
            use_stats: String::default(),
            use_table: String::default(),
            use_trace: String::default(),
            bandwidth: None,
        }
    }
}

impl AlignerBuilder {
    /// Set alignment mode to global (Needleman-Wunsch).
    pub fn global(&mut self) -> &mut Self {
        self.mode = String::from("nw");
        self
    }

    /// Set alignment mode to semi-global
    pub fn semi_global(&mut self) -> &mut Self {
        self.mode = String::from("sg");
        self
    }

    /// Set alignment mode to local (Smith-Waterman).
    pub fn local(&mut self) -> &mut Self {
        self.mode = String::from("sw");
        self
    }

    /// Set solution width (8, 16, 32, or 64 bit). By default, will use sat
    /// (i.e., 8-bit solution width first and falling back to 16-bit if necessary).
    pub fn solution_width(&mut self, solution_width: i32) -> &mut Self {
        self.solution_width = solution_width.to_string();
        self
    }

    /// Set scoring matrix. The default is an identity matrix for DNA sequences.
    /// For more information on creating matrices, see the [Matrix](https://docs.rs/parasail-rs/latest/parasail_rs/struct.Matrix.html) struct.
    /// Default is an identity matrix for DNA sequences.
    pub fn matrix(&mut self, matrix: Matrix) -> &mut Self {
        self.matrix = Arc::new(matrix);
        self
    }

    /// Set gap open penalty.
    /// Note that this should be passed as a positive integer. Default = 5.
    pub fn gap_open(&mut self, gap_open: i32) -> &mut Self {
        self.gap_open = gap_open;
        self
    }

    /// Set gap extend penalty.
    /// Note that this should be passed as a positive integer. Default = 2
    pub fn gap_extend(&mut self, gap_extend: i32) -> &mut Self {
        self.gap_extend = gap_extend;
        self
    }

    /// Set query profile. No query profile is set by default.
    pub fn profile(&mut self, profile: Profile) -> &mut Self {
        self.profile = Arc::new(profile);
        self
    }

    /// Set allowed gaps on query sequence for semi-global alignment.
    /// By default, gaps are allowed at the beginning and end of the query sequence.
    /// Example:
    /// ```rust, no_run
    /// use parasail_rs::Aligner;
    ///
    /// // allow gaps at the beginning of the query sequence
    /// let allow_gaps = vec![String::from("prefix")];
    /// let aligner = Aligner::new().allow_query_gaps(allow_gaps).build();
    /// ```
    pub fn allow_query_gaps(&mut self, allow_gaps: Vec<String>) -> &mut Self {
        self.allow_query_gaps = allow_gaps;
        self
    }

    /// Set allowed gaps on reference sequence for semi-global alignment.
    /// By default, gaps are allowed at the beginning and end of the reference sequence.
    /// Example:
    /// ```rust, no_run
    /// use parasail_rs::Aligner;
    ///
    /// // allow gaps at the beginning of the reference sequence
    /// let allow_gaps = vec![String::from("suffix")];
    /// let aligner = Aligner::new().allow_query_gaps(allow_gaps).build();
    /// ```
    pub fn allow_ref_gaps(&mut self, allow_gaps: Vec<String>) -> &mut Self {
        self.allow_ref_gaps = allow_gaps;
        self
    }

    /// Use striped vectorization method
    pub fn striped(&mut self) -> &mut Self {
        self.vec_strategy = String::from("_striped");
        self
    }

    /// Use scan vectorization method
    pub fn scan(&mut self) -> &mut Self {
        self.vec_strategy = String::from("_scan");
        self
    }

    /// Use diagonal vectorization method
    pub fn diag(&mut self) -> &mut Self {
        self.vec_strategy = String::from("_diag");
        self
    }

    /// Set whether to use statistics for alignment. By default, statistics are
    /// not used. Note that enabling stats and traceback is not supported.
    /// Enabling stats will disable traceback if it is enabled.
    pub fn use_stats(&mut self) -> &mut Self {
        self.use_stats = String::from("_stats");

        // disable traceback if stats are enabled
        if !self.use_trace.is_empty() {
            warn!("Warning: Traceback was enabled previously, but not supported with stats. Disabling traceback");
            self.use_trace = String::default();
        }

        self
    }

    /// Set whether to return the score table. By default, the score table is
    /// not returned. Note that enabling traceback and tables is not supported.
    /// Enabling tables will disable traceback.
    pub fn use_table(&mut self) -> &mut Self {
        self.use_table = String::from("_table");

        // disable traceback if tables are enabled
        if !self.use_trace.is_empty() {
            self.use_trace = String::default();
        }

        self
    }

    /// Set whether to return the last row and column of the score table.
    /// By default, the last row and column are not returned.
    /// Note that if both use_table and use_last_rowcol are set to true, use_table
    /// will be ignored and only the last row and column will be returned.
    pub fn use_last_rowcol(&mut self) -> &mut Self {
        self.use_table = String::from("_rowcol");
        self
    }

    /// Set whether to enable traceback capability. By default, traceback is not enabled.
    /// Note that enabling traceback along with tables or stats is not supported.
    /// Enabling traceback will disable tables and stats if they are enabled.
    pub fn use_trace(&mut self) -> &mut Self {
        self.use_trace = String::from("_trace");

        // disable table if traceback is enabled
        if !self.use_table.is_empty() {
            warn!("Warning: Table was enabled previously, but not supported with traceback. Disabling table");
            self.use_table = String::default();
        }

        // disable stats if traceback is enabled
        if !self.use_stats.is_empty() {
            warn!("Warning: Stats were enabled previously, but not supported with traceback. Disabling stats");
            self.use_stats = String::default();
        }

        self
    }

    /// Helper function for formatting semi-global fn name with correct gap syntax.
    fn get_allowed_gaps(&self, prefix: &str, allowed_gaps_vec: &Vec<String>) -> Vec<String> {
        let mut allowed_gaps = Vec::new();

        if !allowed_gaps_vec.is_empty() {
            if allowed_gaps_vec.contains(&String::from("prefix"))
                && allowed_gaps_vec.contains(&String::from("suffix"))
            {
                allowed_gaps.push(format!("_{prefix}x"));
            } else if allowed_gaps_vec.contains(&String::from("prefix")) {
                allowed_gaps.push(format!("_{prefix}b"));
            } else if allowed_gaps_vec.contains(&String::from("suffix")) {
                allowed_gaps.push(format!("_{prefix}e"));
            }
        }

        allowed_gaps
    }

    /// Get the name of the parasail function to use for alignment.
    fn get_parasail_fn_name(&self) -> CString {
        let mut sg_gaps_fn_part = String::new();
        if self.mode == "sg" {
            let query_gaps_part = self.get_allowed_gaps("q", &self.allow_query_gaps);
            let ref_gaps_part = self.get_allowed_gaps("d", &self.allow_ref_gaps);
            sg_gaps_fn_part = format!("{}{}", query_gaps_part.join(""), ref_gaps_part.join(""));

            if sg_gaps_fn_part == "_qx_dx" {
                sg_gaps_fn_part = String::default();
            }
        }

        let profile: &str;
        let stats: &str;
        if self.profile.is_null() {
            profile = "";
            stats = &self.use_stats;
        } else {
            assert!(
                self.vec_strategy == "_striped" || self.vec_strategy == "_scan",
                "Vectorization strategy must be striped or scan for alignment with a profile."
            );
            profile = "_profile";
            if self.profile.use_stats {
                stats = "_stats";
            } else {
                stats = "";
            }
        }

        CString::new(format!(
            "{}{}{}{}{}{}{}_{}",
            self.mode,
            sg_gaps_fn_part,
            self.use_trace,
            stats,
            self.use_table,
            self.vec_strategy,
            profile,
            self.solution_width
        ))
        .unwrap_or_else(|e| panic!("CString::new failed: {e}"))
    }

    pub fn bandwidth(&mut self, bandwidth: i32) -> &mut Self {
        self.bandwidth = Some(bandwidth);
        self
    }

    /// Build the aligner.
    pub fn build(&mut self) -> Aligner {
        let fn_name = self.get_parasail_fn_name();
        let parasail_fn: AlignerFn;

        if self.profile.is_null() {
            unsafe {
                parasail_fn = AlignerFn::Function(parasail_lookup_function(fn_name.as_ptr()));
            }
        } else {
            unsafe {
                parasail_fn = AlignerFn::PFunction(parasail_lookup_pfunction(fn_name.as_ptr()));
            }
        };

        if parasail_fn.is_none() {
            panic!(
                "Parasail function: {}, not found.",
                fn_name.to_str().unwrap()
            );
        }

        Aligner {
            parasail_fn,
            matrix: Arc::clone(&self.matrix),
            gap_open: self.gap_open,
            gap_extend: self.gap_extend,
            profile: Arc::clone(&self.profile),
            vec_strategy: self.vec_strategy.clone(),
            bandwidth: self.bandwidth,
        }
    }
}

#[derive(Clone)]
/// Aligner struct for sequence alignment
pub struct Aligner {
    parasail_fn: AlignerFn,
    pub matrix: Arc<Matrix>,
    pub gap_open: i32,
    pub gap_extend: i32,
    profile: Arc<Profile>,
    pub vec_strategy: String,
    bandwidth: Option<i32>,
}

impl Aligner {
    /// Create a new default aligner builder.
    pub fn new() -> AlignerBuilder {
        AlignerBuilder::default()
    }

    /// Perform alignment between a query and reference sequence.
    /// If profile was set while building the aligner, pass None as the query
    /// sequence. Otherwise, wrap the query sequence in a Some variant (i.e. Some(query)).
    pub fn align(&self, query: Option<&[u8]>, reference: &[u8]) -> Result<Alignment> {
        let ref_len = reference.len() as i32;
        let reference = CString::new(reference).map_err(Error::NulError)?;

        match self.parasail_fn {
            AlignerFn::Function(f) => {
                assert!(
                    query.is_some(),
                    "Query sequence is required for alignment without a profile."
                );
                let query_raw = query.unwrap();
                let query_len = query_raw.len() as i32;
                let query = CString::new(query_raw).map_err(Error::NulError)?;

                let result = unsafe {
                    // already checked that aligner function f is some variant during build step
                    f.unwrap()(
                        query.as_ptr(),
                        query_len,
                        reference.as_ptr(),
                        ref_len,
                        self.gap_open,
                        self.gap_extend,
                        **self.matrix,
                    )
                };

                Ok(Alignment {
                    inner: result,
                    matrix: **self.matrix,
                    query_len,
                    ref_len,
                })
            }
            AlignerFn::PFunction(f) => {
                // already checked that aligner function f is some variant during build step

                let result = unsafe {
                    f.unwrap()(
                        **self.profile,
                        reference.as_ptr(),
                        ref_len,
                        self.gap_open,
                        self.gap_extend,
                    )
                };

                Ok(Alignment {
                    inner: result,
                    matrix: **self.matrix,
                    query_len: self.profile.query_len,
                    ref_len,
                })
            }
        }
    }

    /// Perform banded global alignment between a query and reference sequence.
    /// Note that this function is not vectorized. However, it may be useful
    /// for aligning large sequences.
    pub fn banded_nw(&self, query: &[u8], reference: &[u8]) -> Result<Alignment> {
        let ref_len = reference.len() as i32;
        let reference = CString::new(reference).map_err(Error::NulError)?;

        let query_len = query.len() as i32;
        let query = CString::new(query).map_err(Error::NulError)?;

        let bandwidth = if let Some(bandwidth) = self.bandwidth {
            bandwidth
        } else {
            return Err(Error::NoBandwidth.into());
        };

        let result = unsafe {
            parasail_nw_banded(
                query.as_ptr(),
                query_len,
                reference.as_ptr(),
                ref_len,
                self.gap_open,
                self.gap_extend,
                bandwidth,
                **self.matrix,
            )
        };

        Ok(Alignment {
            inner: result,
            matrix: **self.matrix,
            query_len,
            ref_len,
        })
    }

    /// Perform Striped Smith-Waterman local alignment using SSE2 instructions.
    pub fn ssw(&self, query: Option<&[u8]>, reference: &[u8]) -> Result<SSWResult> {
        let ref_len = reference.len() as i32;
        let reference = CString::new(reference).map_err(Error::NulError)?;

        let result = if let Some(query_seq) = query {
            let query_len = query_seq.len() as i32;
            let cquery = CString::new(query_seq).map_err(Error::NulError)?;

            unsafe {
                parasail_ssw(
                    cquery.as_ptr(),
                    query_len,
                    reference.as_ptr(),
                    ref_len,
                    self.gap_open,
                    self.gap_extend,
                    self.matrix.inner,
                )
            }
        } else {
            panic!("Query sequence is required for SSW alignment for now.")
            // if self.profile.is_null() {
            //     panic!("Profile is required if no query is provided.")
            // }
            //
            // unsafe {
            //     parasail_ssw_profile(
            //         self.profile.inner,
            //         reference.as_ptr(),
            //         ref_len,
            //         self.gap_open,
            //         self.gap_extend,
            //     )
            // }
        };

        Ok(SSWResult { inner: result })
    }
}

#[doc(hidden)]
unsafe impl Send for Aligner {}
#[doc(hidden)]
unsafe impl Sync for Aligner {}
