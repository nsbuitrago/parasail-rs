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
//! parasail-rs = "0.1"
//! ```
//!
//! ## Examples
//! 
//! ### Basic usage
//! ```rust,no_run
//! use parasail_rs::{ScoringMatrix, Aligner};
//! let matrix = ScoringMatrix::create("ACGT", 1, -1);
//! let vector_strategy = "striped".to_string();
//! let query = b"ACGT";
//! let target = b"ACGTAACGTACA";
//!
//! let aligner = Aligner::new(matrix, 5, 2, vector_strategy);
//! let result = aligner.global(query, target); 
//! ```
//!
//! ### Using query profile
//! ```rust,no_run
//! use parasail_rs::{ScoringMatrix, Aligner, Profile};
//! let matrix = ScoringMatrix::create("ACGT", 1, -1);
//! let vector_strategy = "striped".to_string();
//! let query = b"ACGT";
//! let ref_1 = b"ACGTAACGTACA";
//! let ref_2 = b"TGGCAAGGTAGA";
//!
//! let use_stats = true;
//! let query_profile = Profile::new(query, use_stats, matrix);
//! let aligner = Aligner::with_profile(query_profile, 5, 2, vector_strategy);
//!
//! let result_1 = aligner.global_with_profile(ref_1);
//! let result_2 = aligner.global_with_profile(ref_2);
//!

// use std::marker::PhantomData;

use libparasail_sys::{
    parasail_function_t, parasail_lookup_function, parasail_lookup_pfunction, parasail_matrix_convert_square_to_pssm, parasail_matrix_create, parasail_matrix_free, parasail_matrix_from_file, parasail_matrix_lookup, parasail_matrix_pssm_create, parasail_matrix_t, parasail_pfunction_t, parasail_profile_create_sat, parasail_profile_create_stats_sat, parasail_profile_free, parasail_profile_t, parasail_result_free, parasail_result_get_end_query, parasail_result_get_end_ref, parasail_result_get_length, parasail_result_get_matches, parasail_result_get_score, parasail_result_get_similar, parasail_result_is_banded, parasail_result_is_diag, parasail_result_is_nw, parasail_result_is_saturated, parasail_result_is_scan, parasail_result_is_sg, parasail_result_is_striped, parasail_result_is_sw, parasail_result_t,
};

/// Scoring matrix for sequence alignment.
/// Scores can be created from:
/// - an alphabet and match/mismatch scores
/// - a pre-defined matrix (such as blosum62)
/// - a file containing a scoring matrix (see from_file for format details)
/// - a PSSM (position-specific scoring matrix)
#[derive(Debug, Clone)]
pub struct ScoringMatrix {
    pub scores: *mut parasail_matrix_t,
}

impl ScoringMatrix {
    /// Create a new scoring matrix from an alphabet and match/mismatch scores.
    /// Match score should be a positive integer, while mismatch score should be a negative integer.
    /// The alphabet should be a string containing the characters in the sequences to be aligned.
    pub fn create(alphabet: &str, match_score: i32, mismatch_score: i32) -> ScoringMatrix {
        unsafe {
            let scores = parasail_matrix_create(alphabet.as_ptr() as *const i8, match_score, mismatch_score);
            ScoringMatrix { scores }
        }
    }

    /// Create a new scoring matrix from a pre-defined matrix.
    /// The matrix name should be one of the following:
    /// - blosum{30, 35, 40, 45, 50, 55, 60, 62, 65, 70, 75, 80, 85, 90, 95, 100}
    /// - pam{10-500 in steps of 10}
    pub fn from(matrix_name: &str) -> ScoringMatrix {
        unsafe {
            let scores = parasail_matrix_lookup(matrix_name.as_ptr() as *const i8);
            ScoringMatrix { scores: scores as *mut libparasail_sys::parasail_matrix_t}
        }
    }

    /// Create a new scoring matrix from a file.
    /// Files should contain either square or position-specific scoring matrices.
    /// Examples are direcly from the [Parasail C lib docs](https://github.com/jeffdaily/parasail?tab=readme-ov-file#substitution-matrices).
    /// Square:
    /// ```text
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
    /// ```
    /// PSSM:
    /// ```text
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
    // ```
    pub fn from_file(file: &str) -> ScoringMatrix {
        unsafe {
            let scores = parasail_matrix_from_file(file.as_ptr() as *const i8);
            ScoringMatrix { scores }
        }
    }

    /// Create a new scoring matrix from a PSSM (position-specific scoring matrix).
    pub fn create_pssm(alphabet: &str, values: Vec<i32>, rows: i32) -> ScoringMatrix {
        unsafe {
            let scores = parasail_matrix_pssm_create(alphabet.as_ptr() as *const i8, values.as_ptr(), rows);
            ScoringMatrix { scores }
        }
    }

    /// Convert a square scoring matrix to a PSSM (position-specific scoring matrix).
    pub fn convert_square_to_pssm(&self, pssm_query: &str) -> ScoringMatrix {
        unsafe {
            let scores = parasail_matrix_convert_square_to_pssm(self.scores, pssm_query.as_ptr() as *const i8, pssm_query.len() as i32);
            ScoringMatrix { scores }
        }
    }
}

/// Default scoring matrix is a simple match/mismatch (+1/-1) matrix for DNA sequences.
/// In practice this is never really used, but it's here to allow for default values
/// in the Aligner struct
impl Default for ScoringMatrix {
    fn default() -> Self {
        unsafe {
            let scores = parasail_matrix_create(b"ACGT".as_ptr() as *const i8, 1, -1);
            ScoringMatrix { scores }
        }
    }
}

/// Free the memory used by the scoring matrix when it goes out of scope.
impl Drop for ScoringMatrix {
    fn drop(&mut self) {
        unsafe {
            parasail_matrix_free(self.scores);
        }
    }
}

/// Profile struct to reuse query profiles for alignments. This only applies
/// for striped or scan vector strategies.
#[derive(Debug)]
pub struct Profile {
    pub inner: *mut parasail_profile_t,
}

impl Profile {
    /// Create a new profile from a query sequence, to use with or without stats, and a scoring matrix.
    /// The with_stats should be set to true if you will use an alignment function that returns
    /// statistics. If true, the Profile will use the appropriate parasail functions to allocate
    /// additional data structures required for statistics.
    pub fn new(query: &[u8], with_stats: bool, matrix: ScoringMatrix) -> Profile {
        unsafe {
            let query_len = query.len() as i32;
            match with_stats {
                true => {
                    let query_profile = parasail_profile_create_stats_sat(query.as_ptr() as *const i8, query_len, matrix.scores);
                    Profile { inner: query_profile }
                },
                false => {
                    let query_profile = parasail_profile_create_sat(query.as_ptr() as *const i8, query_len, matrix.scores);
                    Profile { inner: query_profile }
                }
            }
        }
    }
}

/// Free memory used by the profile when it goes out of scope.
impl Drop for Profile {
    fn drop(&mut self) {
        unsafe {
            parasail_profile_free(self.inner);
        }
    }
}

#[derive(Debug)]
pub struct Aligner {
    pub matrix: Option<ScoringMatrix>,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub vec_strategy: String,
    pub query_profile: Option<Profile>,
    pub allow_gaps: Option<Vec<String>>,
}

impl Default for Aligner {
    fn default() -> Self {
        Aligner {
            matrix: Some(ScoringMatrix::default()),
            gap_open: 5,
            gap_extend: 2,
            vec_strategy: String::from("striped"),
            query_profile: None,
            allow_gaps: None,
        }
    }
}

impl Aligner {
    /// Create a new default aligner
    pub fn new() -> Aligner {
        Aligner {
            matrix: Some(ScoringMatrix::default()),
            gap_open: 5,
            gap_extend: 2,
            vec_strategy: String::from("striped"),
            query_profile: None,
            allow_gaps: None,
        }
    }

    /// Set the scoring matrix
    pub fn scoring_matrix(&mut self, matrix: ScoringMatrix) -> &mut Aligner {
        self.matrix = Some(matrix);
        self
    }

    /// Set the gap open penalty for the aligner
    pub fn gap_open_penalty(&mut self, gap_open: i32) -> &mut Aligner {
        self.gap_open = gap_open;
        self
    }

    /// Set the gap extend penalty
    pub fn gap_extend_penalty(&mut self, gap_extend: i32) -> &mut Aligner {
        self.gap_extend = gap_extend;
        self
    }

    /// Set the vector strategy
    pub fn vector_strategy(&mut self, vec_strategy: String) -> &mut Aligner {
        self.vec_strategy = vec_strategy;
        self
    }

    /// Set the query profile
    pub fn query_profile(&mut self, query_profile: Profile) -> &mut Aligner {
        self.query_profile = Some(query_profile);
        self
    }

    /// Set allowed gaps for semi-global alignment. By default, gaps are allwed at the beginning
    /// and end of the query and reference sequences.
    pub fn allow_gaps(&mut self, allow_gaps: Vec<String>) -> &mut Aligner {
        self.allow_gaps = Some(allow_gaps);
        self
    }

    /// Perform alignment between a query and a target sequence. This is a helper method that is
    /// used by the global, semi-global, and local methods. Alternatively, you can call this
    /// directly and pass a mode string of "nw", "sg", or "sw" for global, semi-global, or local, respectively.
    pub fn align(&self, query: Option<&[u8]>, reference: &[u8], mode: String) -> AlignmentResult {
        let reference_len = reference.len() as i32;

        let mut sg_gaps_fn_part = String::new();
        if mode == "sg" {
            if self.allow_gaps.is_some() {
                sg_gaps_fn_part = format!("_{}", self.allow_gaps.as_ref().unwrap().join("_"));
            }
        }

        if self.query_profile.is_some() {
            // use profiles
            // make sure this is only used with striped or scan vector strategies
            assert!(self.vec_strategy == "striped" || self.vec_strategy == "scan");
            let parasail_fn_name = format!("parasail_{}{}_{}_profile_sat", mode, sg_gaps_fn_part, self.vec_strategy);

            unsafe {
                let parasail_fn = parasail_lookup_pfunction(parasail_fn_name.as_ptr() as *const i8);
                if let Some(parasail_fn) = parasail_fn {
                    let result = parasail_fn(self.query_profile.as_ref().unwrap().inner, reference.as_ptr() as *const i8, reference_len, self.gap_open, self.gap_extend);
                    AlignmentResult { result }
                } else {
                    panic!("Invalid alignment method");
                }
            }

        } else {
            // use scoring matrix
            assert!(query.is_some());
            let query = query.unwrap();
            let parasail_fn_name = format!("parasail_{}_{}_profile_sat", mode, self.vec_strategy);

            unsafe {
                let parasail_fn = parasail_lookup_function(parasail_fn_name.as_ptr() as *const i8);
                if let Some(parasail_fn) = parasail_fn {
                    let result = parasail_fn(query.as_ptr() as *const i8, query.len() as i32, reference.as_ptr() as *const i8, reference_len, self.gap_open, self.gap_extend, self.matrix.as_ref().unwrap().scores);
                    AlignmentResult { result }
                } else {
                    panic!("Invalid alignment method");
                }
            }
        }

    }

    pub fn global(&self, query: &[u8], reference: &[u8]) -> AlignmentResult {
        self.align(Some(query), reference, String::from("nw"))
    }

    pub fn global_with_profile(&self, reference: &[u8]) -> AlignmentResult {
        self.align(None, reference, String::from("nw"))
    }

    pub fn semi_global(&self, query: &[u8], reference: &[u8]) -> AlignmentResult {
        self.align(Some(query), reference, String::from("sg"))
    }

    pub fn semi_global_with_profile(&self, reference: &[u8]) -> AlignmentResult {
        self.align(None, reference, String::from("sg"))
    }

    pub fn local(&self, query: &[u8], reference: &[u8]) -> AlignmentResult {
        self.align(Some(query), reference, String::from("sw"))
    }

    pub fn local_with_profile(&self, reference: &[u8]) -> AlignmentResult {
        self.align(None, reference, String::from("sw"))
    }
}

pub enum ParasailFn {
    WithProfile(parasail_pfunction_t),
    WithoutProfile(parasail_function_t),
}

/// Results of a sequence alignment.
#[derive(Debug)]
pub struct AlignmentResult {
    pub result: *mut parasail_result_t,
}

impl AlignmentResult {
    /// Get alignment score.
    pub fn get_score(&self) -> i32 {
        unsafe {
            parasail_result_get_score(self.result)
        }
    }

    pub fn get_end_query(&self) -> i32 {
        unsafe {
            parasail_result_get_end_query(self.result)
        }
    }

    pub fn get_end_ref(&self) -> i32 {
        unsafe {
            parasail_result_get_end_ref(self.result)
        }
    }

    pub fn get_matches(&self) -> i32 {
        unsafe {
            parasail_result_get_matches(self.result)
        }
    }

    pub fn get_similar(&self) -> i32 {
        unsafe {
            parasail_result_get_similar(self.result)
        }
    }

    pub fn get_length(&self) -> i32 {
        unsafe {
            parasail_result_get_length(self.result)
        }
    }

    /// Check if the alignment mode is global.
    pub fn is_global(&self) -> bool {
        unsafe {
            parasail_result_is_nw(self.result) == 1
        }
    }

    /// Check if the alignment mode is semi-global.
    pub fn is_semi_global(&self) -> bool {
        unsafe {
            parasail_result_is_sg(self.result) == 1
        }
    }

    /// Check if the alignment mode is local.
    pub fn is_local(&self) -> bool {
        unsafe {
            parasail_result_is_sw(self.result) == 1
        }
    }

    pub fn is_saturated(&self) -> bool {
        unsafe {
            parasail_result_is_saturated(self.result) == 1
        }
    }

    pub fn is_banded(&self) -> bool {
        unsafe {
            parasail_result_is_banded(self.result) == 1
        }
    }

    /// Check if vector strategy is scan.
    pub fn is_scan(&self) -> bool {
        unsafe {
            parasail_result_is_scan(self.result) == 1
        }
    }

    /// Check if vector strategy is striped.
    pub fn is_striped(&self) -> bool {
        unsafe {
            parasail_result_is_striped(self.result) == 1
        }
    }

    /// Check if vector strategy is diagonal.
    pub fn is_diag(&self) -> bool {
        unsafe {
            parasail_result_is_diag(self.result) == 1
        }
    }
}

/// Free memory used by the alignment result when it goes out of scope.
impl Drop for AlignmentResult {
    fn drop(&mut self) {
        unsafe {
            parasail_result_free(self.result);
        }
    }
}

