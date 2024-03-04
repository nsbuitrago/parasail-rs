//! # Introduction
//!
//! This crate provides safe Rust bindings to the [Parasail](https://github.com/jeffdaily/parasail), a SIMD C library for pairwise sequence alignments.
//!
//! Note that this crate is under development and the bindings may not be complete. For unsafe
//! bindings, see the [libparasail-sys](https://crates.io/crates/libparasail-sys) crate.
//!
//! # Usage
//!
//! Add `parasail-rs` to your `Cargo.toml`:
//!
//! ```toml
//! [dependencies]
//! parasail-rs = "0.2.3"
//! ```
//!
//! ## Examples
//! 
//! ### Basic usage
//! ```rust,no_run
//! use parasail_rs::{Aligner};
//!
//! let query = b"ACGT";
//! let reference = b"ACGT";
//! let aligner = Aligner::new().build();

//! aligner.global(query, reference);
//! ```
//!
//! ### Using query profile
//! ```rust,no_run
//! use parasail_rs::{Matrix, Aligner, Profile};
//!
//! let query = b"ACGT";
//! let ref_1 = b"ACGTAACGTACA";
//! let ref_2 = b"TGGCAAGGTAGA";
//!
//! let use_stats = true;
//! let query_profile = Profile::new(query, use_stats, &Matrix::default());
//! let aligner = Aligner::new()
//!     .profile(query_profile) 
//!     .build();
//!
//! let result_1 = aligner.global_with_profile(ref_1);
//! let result_2 = aligner.global_with_profile(ref_2);
//!

use libparasail_sys::{
    parasail_lookup_function, parasail_lookup_pfunction, parasail_matrix_convert_square_to_pssm, parasail_matrix_create, parasail_matrix_free, parasail_matrix_from_file, parasail_matrix_lookup, parasail_matrix_pssm_create, parasail_matrix_t, parasail_profile_create_sat, parasail_profile_create_stats_sat, parasail_profile_free, parasail_profile_t, parasail_result_free, parasail_result_get_end_query, parasail_result_get_end_ref, parasail_result_get_length, parasail_result_get_length_col, parasail_result_get_length_row, parasail_result_get_length_table, parasail_result_get_matches, parasail_result_get_matches_col, parasail_result_get_matches_row, parasail_result_get_matches_table, parasail_result_get_score, parasail_result_get_score_col, parasail_result_get_score_row, parasail_result_get_score_table, parasail_result_get_similar, parasail_result_get_similar_col, parasail_result_get_similar_row, parasail_result_get_similar_table, parasail_result_get_trace_del_table, parasail_result_get_trace_ins_table, parasail_result_get_trace_table, parasail_result_is_banded, parasail_result_is_blocked, parasail_result_is_diag, parasail_result_is_nw, parasail_result_is_rowcol, parasail_result_is_saturated, parasail_result_is_scan, parasail_result_is_sg, parasail_result_is_stats, parasail_result_is_stats_rowcol, parasail_result_is_stats_table, parasail_result_is_striped, parasail_result_is_sw, parasail_result_is_table, parasail_result_is_trace, parasail_result_t
};

use std::ffi::CString;
use std::ops::Deref;
use std::path::Path;
use std::sync::Arc;
use std::io;

/// Substitution matrix for sequence alignment.
/// Matrices can be created from:
/// - an alphabet and match/mismatch scores
/// - a pre-defined matrix (such as blosum62)
/// - a file containing a substitution matrix (see from_file for format details)
/// - a PSSM (position-specific scoring matrix)
#[derive(Debug, Clone)]
pub struct Matrix {
    inner: *const parasail_matrix_t,
    builtin: bool,
}

impl Matrix {
    /// Create a new scoring matrix from an alphabet and match/mismatch scores.
    /// Note that match score should be a positive integer, while mismatch score should be a negative integer.
    pub fn create(alphabet: &[u8], match_score: i32, mismatch_score: i32) -> Self {
        assert!(match_score >= 0 && mismatch_score <= 0, "Match score should be a positive integer and mismatch score should be a negative integer.");
        unsafe {
            let alphabet = &CString::new(alphabet).unwrap();
            Self { inner: parasail_matrix_create(
                alphabet.as_ptr(),
                match_score,
                mismatch_score),
                builtin: false
            }
        }
    }

    /// Create a new scoring matrix from a pre-defined matrix.
    /// The matrix name should be one of the following:
    /// - blosum{30, 35, 40, 45, 50, 55, 60, 62, 65, 70, 75, 80, 85, 90, 95, 100}
    /// - pam{10-500 in steps of 10}
    pub fn from(matrix_name: &str) -> Self {
        unsafe {
            let matrix = parasail_matrix_lookup(matrix_name.as_ptr() as *const i8);
            Self { inner: matrix, builtin: true }
        }
    }

    /// Create a new scoring matrix from a file.
    /// Files should contain either square or position-specific scoring matrices.
    /// Examples are direcly from the [Parasail C lib docs](https://github.com/jeffdaily/parasail?tab=readme-ov-file#substitution-matrices).
    /// Square:
    ///#
    // # Any line starting with '#' is a comment.
    // #
    // # Needs a row for the alphabet.  First column is a repeat of the
    // # alphabet and assumed to be identical in order to the first alphabet row.
    // #
    // # Last row and column *must* be a non-alphabet character to represent
    // # any input sequence character that is outside of the alphabet.
    // #
    //     A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N   U   *
    // A   5  -4  -4  -4  -4   1   1  -4  -4   1  -4  -1  -1  -1  -2  -4  -5
    // T  -4   5  -4  -4  -4   1  -4   1   1  -4  -1  -4  -1  -1  -2   5  -5
    // G  -4  -4   5  -4   1  -4   1  -4   1  -4  -1  -1  -4  -1  -2  -4  -5
    // C  -4  -4  -4   5   1  -4  -4   1  -4   1  -1  -1  -1  -4  -2  -4  -5
    // S  -4  -4   1   1  -1  -4  -2  -2  -2  -2  -1  -1  -3  -3  -1  -4  -5
    // W   1   1  -4  -4  -4  -1  -2  -2  -2  -2  -3  -3  -1  -1  -1   1  -5
    // R   1  -4   1  -4  -2  -2  -1  -4  -2  -2  -3  -1  -3  -1  -1  -4  -5
    // Y  -4   1  -4   1  -2  -2  -4  -1  -2  -2  -1  -3  -1  -3  -1   1  -5
    // K  -4   1   1  -4  -2  -2  -2  -2  -1  -4  -1  -3  -3  -1  -1   1  -5
    // M   1  -4  -4   1  -2  -2  -2  -2  -4  -1  -3  -1  -1  -3  -1  -4  -5
    // B  -4  -1  -1  -1  -1  -3  -3  -1  -1  -3  -1  -2  -2  -2  -1  -1  -5
    // V  -1  -4  -1  -1  -1  -3  -1  -3  -3  -1  -2  -1  -2  -2  -1  -4  -5
    // H  -1  -1  -4  -1  -3  -1  -3  -1  -3  -1  -2  -2  -1  -2  -1  -1  -5
    // D  -1  -1  -1  -4  -3  -1  -1  -3  -1  -3  -2  -2  -2  -1  -1  -1  -5
    // N  -2  -2  -2  -2  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -2  -5
    // U  -4   5  -4  -4  -4   1  -4   1   1  -4  -1  -4  -1  -1  -2   5  -5
    // *  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5
    /// 
    /// PSSM:
    ///#
    // # Any line starting with '#' is a comment.
    // #
    // # Needs a first row for the alphabet.
    // # First column containing a representative sequence is optional, but included below for an example.
    // #
    //     A   G   I   L   V   M   F   W   P   C   S   T   Y   N   Q   H   K   R   D   E
    // Y  -5  -6   3  -4   2   1   4  -3   0  -5   0  -4   6  0  -5   -4  -5   2  -6  -5
    // S  -1  -5   2  -2   0  -4   1  -6   0   1   3   3  -4  -4  -4  -5  -1   1  -5  -1
    // C  -4  -6  -5  -5  -4  -5  -6  -6  -6  12  -4  -4  -6  -6  -7  -7  -7  -7  -7  -7
    // D  -1  -5  -7  -7  -6  -6  -7  -7  -5  -7  -3  -1  -6   4  -3   3   0  -1   7  -2
    // G   0   4   1  -2  -2  -5  -5  -6  -6   5  -2   0   0   2  -4   3  -5  -5  -5   0
    // C  -4  -6  -5  -5  -4  -5  -6  -6  -6  12  -4  -4  -6  -6  -7  -7  -7  -7  -7  -7
    // L  -4   3  -1   3  -1  -3  -5  -6  -5  -6   0  -4  -5   1   3  -5   1   0  -1  -1
    // K  -2   1   1  -2  -1   3  -5  -6  -5  -5   2   2   0   1   1   1   2  -4  -4   0
    // P  -2   0  -4   0  -2  -4  -5  -5   5  -5  -3  -1   1   1  -3   2  -4  -4   1   3
    // I  -5  -7   7   1   0  -2   3  -5  -6  -5   0  -4  -4  -1  -6   3  -6  -6  -6  -6
    //
    pub fn from_file(file: &str) -> Self {
        let file = CString::new(file).unwrap();
        unsafe {
            Matrix {
                inner: parasail_matrix_from_file(file.as_ptr()), builtin: false
            }
        }
    }

    /// Create a new scoring matrix from a PSSM (position-specific scoring matrix).
    pub fn create_pssm(alphabet: &str, values: Vec<i32>, rows: i32) -> Self {
        let alphabet = CString::new(alphabet).unwrap();
        unsafe {
            Self {
                inner: parasail_matrix_pssm_create(
                    alphabet.as_ptr(),
                    values.as_ptr(),
                    rows
                ),
                builtin: false
            }
        }
    }

    /// Convert a square scoring matrix to a PSSM (position-specific scoring matrix).
    pub fn convert_square_to_pssm(&mut self, pssm_query: &str) {
        let pssm_query_string = CString::new(pssm_query).unwrap();
        unsafe {
            self.inner = parasail_matrix_convert_square_to_pssm(
                    self.inner,
                    pssm_query_string.as_ptr(),
                    pssm_query.len() as i32
                );
        }
    }
}

/// Default scoring matrix is an identity matrix for DNA sequences.
impl Default for Matrix {
    fn default() -> Self {
        Matrix::create(b"ACGT", 1, -1)
    }
}

#[doc(hidden)]
impl Deref for Matrix {
    type Target = *const parasail_matrix_t;
    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

#[doc(hidden)]
impl Drop for Matrix {
    fn drop(&mut self) {

        if self.builtin {
            return;
        }

        unsafe {
            parasail_matrix_free(self.inner as *mut parasail_matrix_t)
        }
    }
}

#[doc(hidden)]
unsafe impl Send for Matrix {}
#[doc(hidden)]
unsafe impl Sync for Matrix {}

// Query profile for sequence alignment
pub struct Profile {
    inner: *mut parasail_profile_t,
    use_stats: bool,
}

impl Profile {
    /// Create a new profile from a query sequence, to use with or without stats, and a scoring matrix.
    /// The with_stats should be set to true if you will use an alignment function that returns
    /// statistics. If true, the Profile will use the appropriate parasail functions to allocate
    /// additional data structures required for statistics.
    pub fn new(query: &[u8], with_stats: bool, matrix: &Matrix) -> Self {
        let query_len = query.len() as i32;
        let query = CString::new(query).unwrap();
        unsafe {
            match with_stats {
                true => {
                    let profile = parasail_profile_create_stats_sat(query.as_ptr(), query_len, **matrix);
                    Profile { inner: profile, use_stats: true }
                },
                false => {
                    let profile = parasail_profile_create_sat(query.as_ptr(), query_len, **matrix);
                    Profile { inner: profile, use_stats: false }
                }
            }
        }
    }
}

/// Default profile is null
impl Default for Profile {
    fn default() -> Self {
        Profile { inner: std::ptr::null_mut(), use_stats: false }
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
            unsafe {
                parasail_profile_free(self.inner)
            }
        }
    }
}

#[doc(hidden)]
unsafe impl Send for Profile {}
#[doc(hidden)]
unsafe impl Sync for Profile {}

// Aligner builder
pub struct AlignerBuilder {
    matrix: Arc<Matrix>,
    gap_open: i32,
    gap_extend: i32,
    profile: Arc<Profile>,
    allow_gaps: Vec<String>,
    vec_strategy: String, 
    use_stats: bool,
    use_table: bool,
    use_last_rowcol: bool,
    use_trace: bool,
}

impl Default for AlignerBuilder {
    fn default() -> Self {
        AlignerBuilder {
            matrix: Matrix::default().into(),
            gap_open: 5,
            gap_extend: 2,
            profile: Profile::default().into(), 
            allow_gaps: Vec::default(),
            vec_strategy: String::from("striped"),
            use_stats: false,
            use_table: false,
            use_last_rowcol: false,
            use_trace: false,
        }
    }
}

impl AlignerBuilder {
    /// Set scoring matrix.
    pub fn matrix(&mut self, matrix: Matrix) -> &mut Self {
        self.matrix = Arc::new(matrix);
        self
    }

    /// Set gap open penalty (Note that this should be passed as a positive integer).
    pub fn gap_open(&mut self, gap_open: i32) -> &mut Self {
        self.gap_open = gap_open;
        self
    }

    /// Set gap extend penalty (Note that this should be passed as a positive integer).
    pub fn gap_extend(&mut self, gap_extend: i32) -> &mut Self {
        self.gap_extend = gap_extend;
        self
    }

    /// Set query profile.
    pub fn profile(&mut self, profile: Profile) -> &mut Self {
        self.profile = Arc::new(profile);
        self
    }

    /// Set allowed gaps at the beginning and end of query/reference sequences 
    /// for semi-global alignment. By default, gaps are allowed at the beginning and end
    /// of both query and reference sequences.
    pub fn allow_gaps(&mut self, allow_gaps: Vec<String>) -> &mut Self {
        self.allow_gaps = allow_gaps;
        self
    }

    /// Set vectorization strategy for alignment. By default, no vectorization
    /// is used.
    pub fn vec_strategy(&mut self, vec_strategy: &str) -> &mut Self {
        self.vec_strategy = vec_strategy.into();
        self
    }

    /// Set whether to use statistics for alignment. By default, statistics are not used.
    pub fn use_stats(&mut self) -> &mut Self {
        self.use_stats = true;
        self
    }

    /// Set whether to return the score table. By default, the score table is not returned.
    pub fn use_table(&mut self) -> &mut Self {
        self.use_table = true;
        self
    }

    /// Set whether to return the last row and column of the score table.
    /// By default, the last row and column are not returned.
    /// Note that if both use_table and use_last_rowcol are set to true, use_table
    /// will be ignored and only the last row and column will be returned.
    pub fn use_last_rowcol(&mut self) -> &mut Self {
        self.use_last_rowcol = true;
        self
    }

    pub fn use_trace(&mut self) -> &mut Self {
        self.use_trace = true;
        self
    }

    /// build aligner
    pub fn build(&mut self) -> Aligner {
        Aligner {
            matrix: Arc::clone(&self.matrix),
            gap_open: self.gap_open,
            gap_extend: self.gap_extend,
            profile: Arc::clone(&self.profile),
            allow_gaps: self.allow_gaps.clone(),
            vec_strategy: self.vec_strategy.clone(),
            use_stats: self.use_stats,
            use_table: self.use_table,
            use_last_rowcol: self.use_last_rowcol,
            use_trace: self.use_trace,
        }
    }
}

/// Aligner struct for sequence alignment
pub struct Aligner {
    matrix: Arc<Matrix>,
    gap_open: i32,
    gap_extend: i32,
    profile: Arc<Profile>,
    allow_gaps: Vec<String>,
    vec_strategy: String,
    use_stats: bool,
    use_table: bool,
    use_last_rowcol: bool,
    use_trace: bool,
}

impl Aligner {
    /// Create a new aligner builder with default fields.
    pub fn new() -> AlignerBuilder {
        AlignerBuilder::default()
    }

    /// Perform alignment between a query and reference sequence.
    /// This is a helper function used by the more specific alignment wrappers.
    /// However, you can call this directly and pass the mode as "nw", "sg", or "sw
    /// for global (Needleman-Wunsch), semi-global, or local (Smith-Watermann) alignment, respectively.
    pub fn align(&self, query: Option<&[u8]>, reference: &[u8], mode: &str) -> AlignResult {
        let ref_len = reference.len() as i32;
        let reference = CString::new(reference).unwrap();
        let mut sg_gaps_fn_part = String::new();
        if mode == "sg" {
            if self.allow_gaps.len() > 0 {
                sg_gaps_fn_part = format!("_{}", self.allow_gaps.join("_"));
            }
        }

        let mut parasail_fn_vec = String::new();
        if self.vec_strategy.len() > 0 {
            parasail_fn_vec = format!("_{}", self.vec_strategy);
        }

        let use_table: &str;
        if self.use_last_rowcol {
            use_table = "_rowcol";
        } else {
            if self.use_table {
                use_table = "_table";
            } else {
                use_table = "";
            }
        }

        let use_trace: &str;
        if self.use_trace {
            use_trace = "_trace";
        } else {
            use_trace = "";
        }

        if self.profile.is_null() {
            // use query
            assert!(query.is_some(), "Query sequence is required for alignment without a profile.");
            let query = query.unwrap();
            let query_len = query.len() as i32;

            let use_stats: &str;
            if self.use_stats {
                use_stats = "_stats";
            } else {
                use_stats = "";
            }

            let parasail_fn_name = CString::new(format!("parasail_{}{}{}{}{}{}_sat", mode, sg_gaps_fn_part, use_trace, use_stats, use_table, parasail_fn_vec))
                .unwrap_or_else(|e| panic!("CString::new failed: {}", e));

            unsafe {
                let parasail_fn = parasail_lookup_function(parasail_fn_name.as_ptr());
                let query = CString::new(query).unwrap();
                if let Some(parasail_fn) = parasail_fn {
                    let result = parasail_fn(
                        query.as_ptr(),
                        query_len,
                        reference.as_ptr(),
                        ref_len,
                        self.gap_open,
                        self.gap_extend,
                        **self.matrix
                    );
                    AlignResult { inner: result }
                } else {
                    panic!("Parasail function: {}, not found.", parasail_fn_name.to_str().unwrap());
                }
            }
        } else {
            // use profile
            
            assert!(self.vec_strategy == "striped" || self.vec_strategy == "scan",
            "Vectorization strategy must be striped or scan for alignment with a profile.");

            let use_stats: &str;
            if self.profile.use_stats {
                use_stats = "_stats";
            } else {
                use_stats = "";
            }

            let parasail_fn_name = CString::new(format!("parasail_{}{}{}{}{}{}_profile_sat", mode, sg_gaps_fn_part, use_trace, use_stats, use_table, parasail_fn_vec))
                .unwrap_or_else(|e| panic!("CString::new failed: {}", e));

            unsafe {
                let parasail_fn = parasail_lookup_pfunction(parasail_fn_name.as_ptr());
                if let Some(parasail_fn) = parasail_fn {
                    let result = parasail_fn(
                        **self.profile,
                        reference.as_ptr(),
                        ref_len,
                        self.gap_open,
                        self.gap_extend,
                    );
                    AlignResult { inner: result }
                } else {
                    panic!("Parasail function: {}, not found.", parasail_fn_name.to_str().unwrap());
                }
            }
        }
    }

    /// Perform global alignment between a query and reference sequence.
    pub fn global(&self, query: &[u8], reference: &[u8]) -> AlignResult {
        self.align(Some(query), reference, "nw")
    }

    /// Perform global alignment using a query profile and reference sequence.
    pub fn global_with_profile(&self, reference: &[u8]) -> AlignResult {
        self.align(None, reference, "nw")
    }

    /// Perform local alignment between a query and reference sequence.
    pub fn local(&self, query: &[u8], reference: &[u8]) -> AlignResult {
        self.align(Some(query), reference, "sw")
    }

    /// Perform local alignment using a query profile and reference sequence.
    pub fn local_with_profile(&self, reference: &[u8]) -> AlignResult {
        self.align(None, reference, "sw")
    }

    /// Perform semi-global alignment between a query and reference sequence.
    pub fn semi_global(&self, query: &[u8], reference: &[u8]) -> AlignResult {
        self.align(Some(query), reference, "sw")
    }

    /// Perform semi-global alignment using a query profile and reference sequence.
    pub fn semi_global_with_profile(&self, reference: &[u8]) -> AlignResult {
        self.align(None, reference, "sw")
    }
}

/// Sequence alignment result.
pub struct AlignResult {
    inner: *mut parasail_result_t
}

impl AlignResult {
    /// Get alignment score.
    pub fn get_score(&self) -> i32 {
        unsafe { parasail_result_get_score(self.inner) }
    }

    /// Get end position of query sequence.
    pub fn get_end_query(&self) -> i32 {
        unsafe {
            parasail_result_get_end_query(self.inner)
        }
    }

    /// Get end position of the reference sequence.
    pub fn get_end_ref(&self) -> i32 {
        unsafe {
            parasail_result_get_end_ref(self.inner)
        }
    }

    /// Get number of matches in the alignment.
    pub fn get_matches(&self) -> Result<i32, io::Error> {
        if self.is_stats() {
            unsafe {
                Ok(parasail_result_get_matches(self.inner))
            }
        } else {
            Err(io::Error::new(io::ErrorKind::Other, "Matches are not available without stats."))
        }
        
    }

    /// 
    pub fn get_similar(&self) -> i32 {
        unsafe {
            parasail_result_get_similar(self.inner)
        }
    }

    /// Get alignment length.
    pub fn get_length(&self) -> Result<i32, io::Error> {
        if self.is_stats() {
            unsafe {
                Ok(parasail_result_get_length(self.inner))
            }
        } else {
            Err(io::Error::new(io::ErrorKind::Other, "Length is not available without stats."))
        }
    }

    /// Get score table.
    pub fn get_score_table(&self) -> Result<i32, io::Error> {
        if self.is_table() || self.is_stats_table() {
            unsafe {
                Ok(*parasail_result_get_score_table(self.inner))
            }
        } else {
            Err(io::Error::new(io::ErrorKind::Other, "Score table is not available without setting use_table"))
        }
    }

    /// Get matches table.
    pub fn get_matches_table(&self) -> Result<i32, io::Error> {
        if self.is_stats_table() {
            unsafe {
                Ok(*parasail_result_get_matches_table(self.inner))
            }
        } else {
            Err(io::Error::new(io::ErrorKind::Other, "Matches table is not available without setting use_stats and use_table"))
        }
    }

    /// Get similar table.
    pub fn get_similar_table(&self) -> Result<i32, io::Error> {
        if self.is_stats_table() {
            unsafe {
                Ok(*parasail_result_get_similar_table(self.inner))
            }
        } else {
            Err(io::Error::new(io::ErrorKind::Other, "Similar table is not available without setting use_stats and use_table"))
        }
    }

    /// Get length table.
    pub fn get_length_table(&self) -> Result<i32, io::Error> {
        if self.is_stats_table() {
            unsafe {
                Ok(*parasail_result_get_length_table(self.inner))
            }
        } else {
            Err(io::Error::new(io::ErrorKind::Other, "Length table is not available without setting use_stats and use_table"))
        }
    }

    /// Get score row.
    pub fn get_score_row(&self) -> Result<i32, io::Error> {
        if self.is_rowcol() || self.is_stats_rowcol() {
            unsafe {
                Ok(*parasail_result_get_score_row(self.inner))
            }
        } else {
            Err(io::Error::new(io::ErrorKind::Other, "Score row is not available without setting use_rowcol"))
        }
    }

    /// Get matches row.
    pub fn get_matches_row(&self) -> Result<i32, io::Error> {
        if self.is_stats_rowcol() {
            unsafe {
                Ok(*parasail_result_get_matches_row(self.inner))
            }
        } else {
            Err(io::Error::new(io::ErrorKind::Other, "Matches row is not available without setting use_stats and use_last_rowcol"))
        }
    }

    /// Get similar row.
    pub fn get_similar_row(&self) -> Result<i32, io::Error> {
        if self.is_stats_rowcol() {
            unsafe {
                Ok(*parasail_result_get_similar_row(self.inner))
            }
        } else {
            Err(io::Error::new(io::ErrorKind::Other, "Similar row is not available without setting use_stats and use_last_rowcol"))
            
        }
    }

    /// Get length row.
    pub fn get_length_row(&self) -> Result<i32, io::Error> {
        if self.is_stats_rowcol() {
            unsafe {
                Ok(*parasail_result_get_length_row(self.inner))
            }
        } else {
            Err(io::Error::new(io::ErrorKind::Other, "Length row is not available without setting use_stats and use_last_rowcol"))
        }
    }

    /// Get score column.
    pub fn get_score_col(&self) -> Result<i32, io::Error> {
        if self.is_rowcol() || self.is_stats_rowcol() {
            unsafe {
                Ok(*parasail_result_get_score_col(self.inner))
            }
        } else {
            Err(io::Error::new(io::ErrorKind::Other, "Score column is not available without setting use_rowcol"))
        }
    }

    /// Get matches column.
    pub fn get_matches_col(&self) -> Result<i32, io::Error> {
        if self.is_stats_rowcol() {
            unsafe {
                Ok(*parasail_result_get_matches_col(self.inner))
            }
        } else {
            Err(io::Error::new(io::ErrorKind::Other, "Matches column is not available without setting use_stats and use_last_rowcol"))
        }
    }

    /// Get similar column.
    pub fn get_similar_col(&self) -> Result<i32, io::Error> {
        if self.is_stats_rowcol() {
            unsafe {
                Ok(*parasail_result_get_similar_col(self.inner))
            }
        } else {
            Err(io::Error::new(io::ErrorKind::Other, "Similar column is not available without setting use_stats and use_last_rowcol"))
        }
    }

    /// Get length column
    pub fn get_length_col(&self) -> Result<i32, io::Error> {
        if self.is_stats_rowcol() {
            unsafe {
                Ok(*parasail_result_get_length_col(self.inner))
            }
        } else {
            Err(io::Error::new(io::ErrorKind::Other, "Length column is not available without setting use_stats and use_last_rowcol"))
        }
    }

    /// Get trace table.
    pub fn get_trace_table(&self) -> Result<i32, io::Error> {
        if self.is_trace() {
            unsafe {
                Ok(*parasail_result_get_trace_table(self.inner))
            }
        } else {
            Err(io::Error::new(io::ErrorKind::Other, "Trace table is not available without setting use_trace"))
        }
    }

    /// Get trace insertion table.
    pub fn get_trace_ins_table(&self) -> Result<i32, io::Error> {
        if self.is_trace() {
            unsafe {
                Ok(*parasail_result_get_trace_ins_table(self.inner))
            }
        } else {
            Err(io::Error::new(io::ErrorKind::Other, "Trace insertion table is not available without setting use_trace"))
        }
    }

    /// Get trace deletion table.
    pub fn get_trace_del_table(&self) -> Result<i32, io::Error> {
        if self.is_trace() {
            unsafe {
                Ok(*parasail_result_get_trace_del_table(self.inner))
            }
        } else {
            Err(io::Error::new(io::ErrorKind::Other, "Trace insertion table is not available without setting use_trace"))
        }
    }

    /// Check if the alignment mode is global.
    pub fn is_global(&self) -> bool {
        unsafe {
            parasail_result_is_nw(self.inner) != 0
        }
    }

    /// Check if the alignment mode is semi-global.
    pub fn is_semi_global(&self) -> bool {
        unsafe {
            parasail_result_is_sg(self.inner) != 0
        }
    }

    /// Check if the alignment mode is local.
    pub fn is_local(&self) -> bool {
        
        unsafe {
            parasail_result_is_sw(self.inner) != 0
        }
    }

    /// Check if the solution width is saturated (i.e., using 8-bit solution width first and
    /// falling back to 16-bit if necessary).
    pub fn is_saturated(&self) -> bool {
        unsafe {
            parasail_result_is_saturated(self.inner) != 0
        }
    }

    /// Check if banded alignment is used. 
    pub fn is_banded(&self) -> bool {
        unsafe {
            parasail_result_is_banded(self.inner) != 0
        }
    }

    /// Check if vector strategy is scan.
    pub fn is_scan(&self) -> bool {
        unsafe {
            parasail_result_is_scan(self.inner) != 0
        }
    }

    /// Check if vector strategy is striped.
    pub fn is_striped(&self) -> bool {
        unsafe {
            parasail_result_is_striped(self.inner) != 0
        }
    }

    /// Check if vector strategy is diagonal.
    pub fn is_diag(&self) -> bool {
        unsafe {
            parasail_result_is_diag(self.inner) != 0
        }
    }

    pub fn is_blocked(&self) -> bool {
        unsafe {
            parasail_result_is_blocked(self.inner) != 0
        }
    }

    /// Check if statistics are returned from alignment.
    pub fn is_stats(&self) -> bool {
        unsafe {
            parasail_result_is_stats(self.inner) != 0
        }
    }

    /// Check if result is a stats table
    pub fn is_stats_table(&self) -> bool {
        unsafe {
            parasail_result_is_stats_table(self.inner) != 0
        }
    }

    /// Check if result is a table
    pub fn is_table(&self) -> bool {
        unsafe {
            parasail_result_is_table(self.inner) != 0
        }
    }

    /// Check if result is a last row and column of table
    pub fn is_rowcol(&self) -> bool {
        unsafe {
            parasail_result_is_rowcol(self.inner) != 0
        }
    }

    /// Check if result is a row and column of table with additional statistics.
    pub fn is_stats_rowcol(&self) -> bool {
        unsafe {
            parasail_result_is_stats_rowcol(self.inner) != 0
        }
    }

    /// Check if is trace.
    pub fn is_trace(&self) -> bool {
        unsafe {
            parasail_result_is_trace(self.inner) != 0
        }
    }
}

#[doc(hidden)]
impl Drop for AlignResult {
    fn drop(&mut self) {
        unsafe {
            parasail_result_free(self.inner);
        }
    }
}

