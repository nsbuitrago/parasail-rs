//! # Introduction
//!
//! This crate provides safe Rust bindings to
//! [Parasail](https://github.com/jeffdaily/parasail), a SIMD pairwise sequence
//! alignment library.
//!
//! For unsafe bindings, see the
//! [libparasail-sys](https://crates.io/crates/libparasail-sys) crate.
//!
//! # Usage
//!
//! ## Examples
//!
//! ### Basic usage
//! ```rust,no_run
//! use parasail_rs::{Aligner};
//!
//! fn align() -> Result<(), Box<dyn std::error::Error>> {
//!     let query = b"ACGT";
//!     let reference = b"ACGT";
//!     let aligner = Aligner::new().build();

//!     let result = aligner.align(Some(query), reference)?;
//!     println!("Alignment Score: {}", result.get_score());
//!
//!     Ok(())
//! }
//! ```
//!
//! ### Using query profile
//!
//! When using striped or scan vectorization strategies, some performance may
//! be gained by reusing the query sequence. This can be done by creating a
//! query profile and reusing it for multiple alignments.
//!
//! ```rust,no_run
//! use parasail_rs::{Matrix, Aligner, Profile};
//!
//! fn align() -> Result<(), Box<dyn std::error::Error>> {
//!
//!     let query = b"ACGT";
//!     let ref_1 = b"ACGTAACGTACA";
//!     let ref_2 = b"TGGCAAGGTAGA";
//!
//!     let use_stats = true;
//!     let query_profile = Profile::new(query, use_stats, &Matrix::default())?;
//!     let aligner = Aligner::new()
//!         .profile(query_profile)
//!         .build();
//!
//!     let result_1 = aligner.align(None, ref_1)?;
//!     let result_2 = aligner.align(None, ref_2)?;
//!
//!     println!("Score 1: {}", result_1.get_score());
//!     println!("Score 2: {}", result_2.get_score());
//!     Ok(())
//! }
//!

use libc::c_char;
use libparasail_sys::{
    parasail_cigar_decode, parasail_cigar_free, parasail_cigar_t, parasail_lookup_function,
    parasail_lookup_pfunction, parasail_matrix_convert_square_to_pssm, parasail_matrix_copy,
    parasail_matrix_create, parasail_matrix_free, parasail_matrix_from_file,
    parasail_matrix_lookup, parasail_matrix_pssm_create, parasail_matrix_set_value,
    parasail_matrix_t, parasail_nw_banded, parasail_profile_create_sat,
    parasail_profile_create_stats_sat, parasail_profile_free, parasail_profile_t,
    parasail_result_free, parasail_result_get_cigar, parasail_result_get_end_query,
    parasail_result_get_end_ref, parasail_result_get_length, parasail_result_get_length_col,
    parasail_result_get_length_row, parasail_result_get_length_table, parasail_result_get_matches,
    parasail_result_get_matches_col, parasail_result_get_matches_row,
    parasail_result_get_matches_table, parasail_result_get_score, parasail_result_get_score_col,
    parasail_result_get_score_row, parasail_result_get_score_table, parasail_result_get_similar,
    parasail_result_get_similar_col, parasail_result_get_similar_row,
    parasail_result_get_similar_table, parasail_result_get_trace_table,
    parasail_result_get_traceback, parasail_result_is_banded, parasail_result_is_blocked,
    parasail_result_is_diag, parasail_result_is_nw, parasail_result_is_rowcol,
    parasail_result_is_saturated, parasail_result_is_scan, parasail_result_is_sg,
    parasail_result_is_stats, parasail_result_is_stats_rowcol, parasail_result_is_stats_table,
    parasail_result_is_striped, parasail_result_is_sw, parasail_result_is_table,
    parasail_result_is_trace, parasail_result_ssw_free, parasail_result_ssw_t, parasail_result_t,
    parasail_ssw, parasail_traceback_generic,
};

use log::warn;
use std::ffi::{CString, IntoStringError, NulError};
use std::fmt::Display;
use std::ops::Deref;
use std::path::Path;
use std::slice;
use std::sync::Arc;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum MatrixError {
    #[error("Error creating matrix: {0}")]
    CreateErr(#[from] NulError),
    #[error("Error creating matrix from lookup: {0} matrix not found.")]
    LookupErr(String),
    #[error("Error creating matrix from file: {0}")]
    FromFileErr(String),
    #[error("File not found: {0}")]
    FileNotFound(String),
    #[error("Error creating PSSM matrix. Invalid alphabet, values, or rows.")]
    CreatePssmErr,
    #[error("Error creating matrix. Matrix has not been created yet. Consider using the `create` method.")]
    NullMatrix,
    #[error("Error converting matrix to PSSM. Matrix is not square.")]
    NotSquare,
    #[error("Error setting value on substitution matrix: Matrix must be a user matrix and not builtin. Consider using the `create` method")]
    NotBuiltIn,
    #[error("Error setting value on substitution matrix: Index ({0},{1}) out of range ((0,0),({2},{2}))")]
    InvalidIndex(i32, i32, i32),
}

#[derive(Error, Debug)]
pub enum ProfileError {
    #[error("Error creating profile: {0}")]
    CreateErr(#[from] NulError),
    #[error("Error creating profile. Null profile returned from parasail.")]
    NullProfile,
}

#[derive(Error, Debug)]
pub enum AlignError {
    #[error("Alignment initialization error: {0}")]
    AlignInitErr(#[from] NulError),
    #[error("No bandwith set for banded alignment.")]
    NoBandwith,
}

#[derive(Error, Debug)]
pub enum AlignResultError {
    #[error("Error getting result from {0}. Stats must be initialized. Consider using `use_stats` method on AlignerBuilder.")]
    NoStats(String),
    #[error("Error getting result from {0}. Table must be enabled. Consider using `use_table` method on AlignerBuilder.")]
    NoTable(String),
    #[error("Error getting result from {0}. Table and stats must be enabled. Consider chaining `use_stats` and use_table methods on AlignerBuilder.")]
    NoStatsTable(String),
    #[error("Error getting result from {0}. Last row and col must be enabled. Consider using `use_rowcol` method on AlignerBuilder.")]
    NoRowCol(String),
    #[error("Error getting result from {0}. Traceback must be enabled. Consider using `use_trace` method on AlignerBuilder.")]
    NoTrace(String),
    #[error("Error converting CIGAR string to Rust string: {0}")]
    CigarToStringErr(#[from] IntoStringError),
    #[error("Error creating new CString: {0}")]
    NewCStringErr(#[from] NulError),
    #[error("No bandwith set for banded alignment.")]
    NoBandwith,
}

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
    /// Note that match score should be a positive integer, while mismatch score
    /// should be a negative integer.
    pub fn create(
        alphabet: &[u8],
        match_score: i32,
        mismatch_score: i32,
    ) -> Result<Self, MatrixError> {
        assert!(match_score >= 0 && mismatch_score <= 0, "Match score should be a positive integer and mismatch score should be a negative integer.");
        assert!(alphabet.len() > 0, "Alphabet should not be empty.");
        unsafe {
            let alphabet = &CString::new(alphabet)?;
            Ok(Self {
                inner: parasail_matrix_create(alphabet.as_ptr(), match_score, mismatch_score),
                builtin: false,
            })
        }
    }

    /// Create a new scoring matrix from a pre-defined matrix.
    /// The matrix name should be one of the following:
    /// - blosum{30, 35, 40, 45, 50, 55, 60, 62, 65, 70, 75, 80, 85, 90, 95, 100}
    /// - pam{10-500} (in steps of 10, i.e., pam10, pam20, ... pam500).
    pub fn from(matrix_name: &str) -> Result<Self, MatrixError> {
        assert!(!matrix_name.is_empty(), "Matrix name should not be empty.");
        let matrix: *const parasail_matrix_t;
        unsafe {
            let matrix_name = CString::new(matrix_name)?;
            matrix = parasail_matrix_lookup(matrix_name.as_ptr());
        }

        if matrix.is_null() {
            return Err(MatrixError::LookupErr(matrix_name.to_string()));
        }

        Ok(Self {
            inner: matrix,
            builtin: true,
        })
    }

    /// Create a new scoring matrix from a file.
    /// Files should contain either square or position-specific scoring matrices.
    /// Examples are direcly from the [Parasail C lib docs](https://github.com/jeffdaily/parasail?tab=readme-ov-file#substitution-matrices).
    /// Square:
    /// ```plaintext
    ///#
    /// # Any line starting with '#' is a comment.
    /// #
    /// # Needs a row for the alphabet.  First column is a repeat of the
    /// # alphabet and assumed to be identical in order to the first alphabet row.
    /// #
    /// # Last row and column *must* be a non-alphabet character to represent
    /// # any input sequence character that is outside of the alphabet.
    /// #
    ///     A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N   U   *
    /// A   5  -4  -4  -4  -4   1   1  -4  -4   1  -4  -1  -1  -1  -2  -4  -5
    /// T  -4   5  -4  -4  -4   1  -4   1   1  -4  -1  -4  -1  -1  -2   5  -5
    /// G  -4  -4   5  -4   1  -4   1  -4   1  -4  -1  -1  -4  -1  -2  -4  -5
    /// C  -4  -4  -4   5   1  -4  -4   1  -4   1  -1  -1  -1  -4  -2  -4  -5
    /// S  -4  -4   1   1  -1  -4  -2  -2  -2  -2  -1  -1  -3  -3  -1  -4  -5
    /// W   1   1  -4  -4  -4  -1  -2  -2  -2  -2  -3  -3  -1  -1  -1   1  -5
    /// R   1  -4   1  -4  -2  -2  -1  -4  -2  -2  -3  -1  -3  -1  -1  -4  -5
    /// Y  -4   1  -4   1  -2  -2  -4  -1  -2  -2  -1  -3  -1  -3  -1   1  -5
    /// K  -4   1   1  -4  -2  -2  -2  -2  -1  -4  -1  -3  -3  -1  -1   1  -5
    /// M   1  -4  -4   1  -2  -2  -2  -2  -4  -1  -3  -1  -1  -3  -1  -4  -5
    /// B  -4  -1  -1  -1  -1  -3  -3  -1  -1  -3  -1  -2  -2  -2  -1  -1  -5
    /// V  -1  -4  -1  -1  -1  -3  -1  -3  -3  -1  -2  -1  -2  -2  -1  -4  -5
    /// H  -1  -1  -4  -1  -3  -1  -3  -1  -3  -1  -2  -2  -1  -2  -1  -1  -5
    /// D  -1  -1  -1  -4  -3  -1  -1  -3  -1  -3  -2  -2  -2  -1  -1  -1  -5
    /// N  -2  -2  -2  -2  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -2  -5
    /// U  -4   5  -4  -4  -4   1  -4   1   1  -4  -1  -4  -1  -1  -2   5  -5
    /// *  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5
    /// ```
    ///
    /// PSSM:
    /// ```plaintext
    ///#
    /// # Any line starting with '#' is a comment.
    /// #
    /// # Needs a first row for the alphabet.
    /// # First column containing a representative sequence is optional, but included below for an example.
    /// #
    ///     A   G   I   L   V   M   F   W   P   C   S   T   Y   N   Q   H   K   R   D   E
    /// Y  -5  -6   3  -4   2   1   4  -3   0  -5   0  -4   6  0  -5   -4  -5   2  -6  -5
    /// S  -1  -5   2  -2   0  -4   1  -6   0   1   3   3  -4  -4  -4  -5  -1   1  -5  -1
    /// C  -4  -6  -5  -5  -4  -5  -6  -6  -6  12  -4  -4  -6  -6  -7  -7  -7  -7  -7  -7
    /// D  -1  -5  -7  -7  -6  -6  -7  -7  -5  -7  -3  -1  -6   4  -3   3   0  -1   7  -2
    /// G   0   4   1  -2  -2  -5  -5  -6  -6   5  -2   0   0   2  -4   3  -5  -5  -5   0
    /// C  -4  -6  -5  -5  -4  -5  -6  -6  -6  12  -4  -4  -6  -6  -7  -7  -7  -7  -7  -7
    /// L  -4   3  -1   3  -1  -3  -5  -6  -5  -6   0  -4  -5   1   3  -5   1   0  -1  -1
    /// K  -2   1   1  -2  -1   3  -5  -6  -5  -5   2   2   0   1   1   1   2  -4  -4   0
    /// P  -2   0  -4   0  -2  -4  -5  -5   5  -5  -3  -1   1   1  -3   2  -4  -4   1   3
    /// I  -5  -7   7   1   0  -2   3  -5  -6  -5   0  -4  -4  -1  -6   3  -6  -6  -6  -6
    /// ```
    pub fn from_file(file: &str) -> Result<Self, MatrixError> {
        let filepath = Path::new(file);
        if !filepath.exists() {
            return Err(MatrixError::FileNotFound(
                filepath.to_str().unwrap().to_string(),
            ));
        }

        let file = CString::new(file)?;

        unsafe {
            let matrix = parasail_matrix_from_file(file.as_ptr());

            if matrix.is_null() {
                return Err(MatrixError::FromFileErr(
                    filepath.to_str().unwrap().to_string(),
                ));
            }

            Ok(Self {
                inner: parasail_matrix_from_file(file.as_ptr()),
                builtin: false,
            })
        }
    }

    /// Create a new scoring matrix from a PSSM (position-specific scoring matrix).
    pub fn create_pssm(alphabet: &str, values: Vec<i32>, rows: i32) -> Result<Self, MatrixError> {
        let alphabet = CString::new(alphabet)?;

        unsafe {
            let matrix = parasail_matrix_pssm_create(alphabet.as_ptr(), values.as_ptr(), rows);

            if matrix.is_null() {
                return Err(MatrixError::CreatePssmErr);
            }

            Ok(Self {
                inner: matrix,
                builtin: false,
            })
        }
    }

    /// Convert a square scoring matrix to a PSSM (position-specific scoring matrix).
    pub fn to_pssm(self, pssm_query: &[u8]) -> Result<Matrix, MatrixError> {
        assert!(
            !pssm_query.is_empty(),
            "PSSM query sequence should not be empty."
        );
        let pssm_query_string = CString::new(pssm_query)?;

        unsafe {
            let matrix = parasail_matrix_copy(self.inner);
            if matrix.is_null() {
                return Err(MatrixError::NullMatrix);
            }

            if (*self.inner).type_ != 0 {
                return Err(MatrixError::NotSquare);
            }

            let converted_matrix = parasail_matrix_convert_square_to_pssm(
                matrix,
                pssm_query_string.as_ptr(),
                pssm_query.len() as i32,
            );

            if converted_matrix.is_null() {
                panic!("Erro converting matrix to PSSM. Invalid query sequence.")
            }

            Ok(Matrix {
                inner: converted_matrix,
                builtin: self.builtin,
            })
        }
    }

    /// Set value at a given row and column in a user defined substitution matrix.
    pub fn set_value(&mut self, row: i32, col: i32, value: i32) -> Result<(), MatrixError> {
        if self.builtin {
            return Err(MatrixError::NotBuiltIn);
        }

        unsafe {
            let size = (*self.inner).size - 2;

            if size < 0 {
                return Err(MatrixError::NullMatrix);
            }

            if row < 0 || row > size || col < 0 || col > size {
                return Err(MatrixError::InvalidIndex(row, col, size));
            }

            parasail_matrix_set_value(self.inner.cast_mut(), row, col, value);
        }

        Ok(())
    }
}

/// Default scoring matrix is an identity matrix for DNA sequences.
impl Default for Matrix {
    fn default() -> Self {
        Matrix::create(b"ACGTA", 1, -1).unwrap()
    }
}

#[doc(hidden)]
impl Display for Matrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        unsafe {
            let length = (*self.inner).length as usize;
            let size = (*self.inner).size as usize;
            let matrix = slice::from_raw_parts((*self.inner).matrix, length * size);
            for i in 0..length {
                for j in 0..size {
                    write!(f, "{} ", matrix[i * size + j])?;
                }
                writeln!(f)?;
            }
        }
        Ok(())
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

        unsafe { parasail_matrix_free(self.inner as *mut parasail_matrix_t) }
    }
}

#[doc(hidden)]
unsafe impl Send for Matrix {}
#[doc(hidden)]
unsafe impl Sync for Matrix {}

/// Query profile for sequence alignment
pub struct Profile {
    inner: *mut parasail_profile_t,
    use_stats: bool,
}

impl Profile {
    /// Create a new profile from a query sequence, to use with or without stats, and a scoring matrix.
    /// The with_stats should be set to true if you will use an alignment function that returns
    /// statistics. If true, the Profile will use the appropriate parasail functions to allocate
    /// additional data structures required for statistics.
    pub fn new(query: &[u8], with_stats: bool, matrix: &Matrix) -> Result<Self, ProfileError> {
        let query_len = query.len() as i32;
        if query_len == 0 {
            panic!("Query sequence is empty.");
        }
        let query = CString::new(query)?;

        unsafe {
            match with_stats {
                true => {
                    let profile =
                        parasail_profile_create_stats_sat(query.as_ptr(), query_len, **matrix);
                    if profile.is_null() {
                        return Err(ProfileError::NullProfile);
                    }

                    Ok(Profile {
                        inner: profile,
                        use_stats: true,
                    })
                }
                false => {
                    let profile = parasail_profile_create_sat(query.as_ptr(), query_len, **matrix);
                    if profile.is_null() {
                        return Err(ProfileError::NullProfile);
                    }

                    Ok(Profile {
                        inner: profile,
                        use_stats: false,
                    })
                }
            }
        }
    }

    //pub fn ssw_init(query: &[u8], matrix: &Matrix, score_size: i8) -> Result<Self, ProfileError> {
    //    let query_len = query.len() as i32;
    //    if query_len == 0 {
    //        panic!("Query sequence is empty.");
    //    }
    //    let query = CString::new(query)?;
    //
    //    unsafe {
    //        let profile = parasail_ssw_init(query.as_ptr(), query_len, **matrix, score_size);
    //        if profile.is_null() {
    //            return Err(ProfileError::NullProfile);
    //        }
    //
    //        Ok(Profile {
    //            inner: profile,
    //            use_stats: false,
    //        })
    //    }
    //}
}

/// Default profile is a null pointer
// This is for cases where the profile is not used by the Aligner but we need
// some default anyway. We could probably also not have this default and
// just wrap the Profile in an Option
impl Default for Profile {
    fn default() -> Self {
        Profile {
            inner: std::ptr::null_mut(),
            use_stats: false,
        }
    }
}

#[doc(hidden)]
impl Deref for Profile {
    type Target = *mut parasail_profile_t;
    fn deref(&self) -> &Self::Target {
        // also check if inner is null
        // and otherwise just return nothing
        if !self.inner.is_null() {
            return;
        }

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

/// Parasail alignment function type.
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
    bandwith: Option<i32>,
}

/// Default aligner uses global alignment with an identity matrix for DNA
/// sequences and no gap penalties. No profile, trace, table, or stats options
/// are set. Vectorization strategy is set to striped by default.
impl Default for AlignerBuilder {
    fn default() -> Self {
        AlignerBuilder {
            mode: String::from("nw"),
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
            bandwith: None,
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

    /// Set alignment mode to local (Smith-Watermann).
    pub fn local(&mut self) -> &mut Self {
        self.mode = String::from("sw");
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

        if allowed_gaps_vec.len() > 0 {
            if allowed_gaps_vec.contains(&String::from("prefix"))
                && allowed_gaps_vec.contains(&String::from("suffix"))
            {
                allowed_gaps.push(format!("_{}x", prefix));
            } else if allowed_gaps_vec.contains(&String::from("prefix")) {
                allowed_gaps.push(format!("_{}b", prefix));
            } else if allowed_gaps_vec.contains(&String::from("suffix")) {
                allowed_gaps.push(format!("_{}e", prefix));
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

        let parasail_fn_name = CString::new(format!(
            "{}{}{}{}{}{}{}_sat",
            self.mode,
            sg_gaps_fn_part,
            self.use_trace,
            stats,
            self.use_table,
            self.vec_strategy,
            profile
        ))
        .unwrap_or_else(|e| panic!("CString::new failed: {}", e));

        parasail_fn_name
    }

    pub fn bandwith(&mut self, bandwith: i32) -> &mut Self {
        self.bandwith = Some(bandwith);
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
            bandwith: self.bandwith,
        }
    }
}

/// Aligner struct for sequence alignment
pub struct Aligner {
    parasail_fn: AlignerFn,
    pub matrix: Arc<Matrix>,
    pub gap_open: i32,
    pub gap_extend: i32,
    profile: Arc<Profile>,
    pub vec_strategy: String,
    bandwith: Option<i32>,
}

impl Aligner {
    /// Create a new default aligner builder.
    pub fn new() -> AlignerBuilder {
        AlignerBuilder::default()
    }

    /// Perform alignment between a query and reference sequence.
    /// If profile was set while building the aligner, pass None as the query
    /// sequence. Otherwise, wrap the query sequence in a Some variant (i.e. Some(query)).
    pub fn align(&self, query: Option<&[u8]>, reference: &[u8]) -> Result<AlignResult, AlignError> {
        let ref_len = reference.len() as i32;
        let reference = CString::new(reference)?;

        match self.parasail_fn {
            AlignerFn::Function(f) => {
                assert!(
                    query.is_some(),
                    "Query sequence is required for alignment without a profile."
                );
                let query = query.unwrap();
                let query_len = query.len() as i32;

                unsafe {
                    let query = CString::new(query)?;
                    // already checked that aligner function f is some variant during build step
                    let result = f.unwrap()(
                        query.as_ptr(),
                        query_len,
                        reference.as_ptr(),
                        ref_len,
                        self.gap_open,
                        self.gap_extend,
                        **self.matrix,
                    );

                    Ok(AlignResult {
                        inner: result,
                        matrix: **self.matrix,
                    })
                }
            }
            AlignerFn::PFunction(f) => unsafe {
                // already checked that aligner function f is some variant during build step
                let result = f.unwrap()(
                    **self.profile,
                    reference.as_ptr(),
                    ref_len,
                    self.gap_open,
                    self.gap_extend,
                );

                Ok(AlignResult {
                    inner: result,
                    matrix: **self.matrix,
                })
            },
        }
    }

    /// Peform banded global alignment between a query and reference sequence.
    /// Note that this function is not vectorized. However, it may be useful
    /// for aligning large sequences.
    pub fn banded_nw(&self, query: &[u8], reference: &[u8]) -> Result<AlignResult, AlignError> {
        let ref_len = reference.len() as i32;
        let reference = CString::new(reference)?;

        let query_len = query.len() as i32;
        let query = CString::new(query)?;

        let bandwith = if let Some(bandwith) = self.bandwith {
            bandwith
        } else {
            return Err(AlignError::NoBandwith);
        };

        unsafe {
            let result = parasail_nw_banded(
                query.as_ptr(),
                query_len,
                reference.as_ptr(),
                ref_len,
                self.gap_open,
                self.gap_extend,
                bandwith,
                **self.matrix,
            );

            Ok(AlignResult {
                inner: result,
                matrix: **self.matrix,
            })
        }
    }

    /// Perform Striped Smith-Waterman local alignment using SSE2 instructions.
    pub fn ssw(&self, query: Option<&[u8]>, reference: &[u8]) -> Result<SSWResult, AlignError> {
        let ref_len = reference.len() as i32;
        let reference = CString::new(reference)?;

        if query.is_some() {
            let query = query.unwrap();
            let query_len = query.len() as i32;

            unsafe {
                let query = CString::new(query)?;
                let result = parasail_ssw(
                    query.as_ptr(),
                    query_len,
                    reference.as_ptr(),
                    ref_len,
                    self.gap_open,
                    self.gap_extend,
                    **self.matrix,
                );

                Ok(SSWResult {
                    inner: result, // matrix: **self.matrix,
                })
            }
        } else {
            panic!("Query sequence is required for SSW alignment for now.");
            //unsafe {
            //    // need to check if the profile is of type parasail_profile_t
            //    let result = parasail_ssw_profile(
            //        **self.profile,
            //        reference.as_ptr(),
            //        ref_len,
            //        self.gap_open,
            //        self.gap_extend,
            //    );
            //
            //    Ok(SSWResult {
            //        inner: result, //matrix: **self.matrix,
            //    })
            //}
        }
    }
}

#[doc(hidden)]
unsafe impl Send for Aligner {}
#[doc(hidden)]
unsafe impl Sync for Aligner {}

/// CIGAR string for sequence alignment.
struct CigarString {
    inner: *mut parasail_cigar_t,
}

#[doc(hidden)]
impl Drop for CigarString {
    fn drop(&mut self) {
        unsafe {
            parasail_cigar_free(self.inner);
        }
    }
}

/// Traceback for sequence alignment.
pub struct Traceback {
    pub query: String,
    pub comparison: String,
    pub reference: String,
}

/// Sequence alignment result.
#[derive(Debug, Clone)]
pub struct AlignResult {
    inner: *mut parasail_result_t,
    matrix: *const parasail_matrix_t,
}

impl AlignResult {
    /// Get alignment score.
    pub fn get_score(&self) -> i32 {
        unsafe { parasail_result_get_score(self.inner) }
    }

    /// Get end position of query sequence.
    pub fn get_end_query(&self) -> i32 {
        unsafe { parasail_result_get_end_query(self.inner) }
    }

    /// Get end position of the reference sequence.
    pub fn get_end_ref(&self) -> i32 {
        unsafe { parasail_result_get_end_ref(self.inner) }
    }

    /// Get number of matches in the alignment.
    pub fn get_matches(&self) -> Result<i32, AlignResultError> {
        if self.is_stats() {
            unsafe { Ok(parasail_result_get_matches(self.inner)) }
        } else {
            Err(AlignResultError::NoStats(String::from("get_matches()")))
        }
    }

    pub fn get_similar(&self) -> i32 {
        unsafe { parasail_result_get_similar(self.inner) }
    }

    /// Get alignment length.
    pub fn get_length(&self) -> Result<i32, AlignResultError> {
        if self.is_stats() {
            unsafe { Ok(parasail_result_get_length(self.inner)) }
        } else {
            Err(AlignResultError::NoStats(String::from("get_length()")))
        }
    }

    /// Get score table.
    pub fn get_score_table(&self) -> Result<i32, AlignResultError> {
        if self.is_table() || self.is_stats_table() {
            unsafe { Ok(*parasail_result_get_score_table(self.inner)) }
        } else {
            Err(AlignResultError::NoTable(String::from("get_score_table()")))
        }
    }

    /// Get matches table.
    pub fn get_matches_table(&self) -> Result<i32, AlignResultError> {
        if self.is_stats_table() {
            unsafe { Ok(*parasail_result_get_matches_table(self.inner)) }
        } else {
            Err(AlignResultError::NoStatsTable(String::from(
                "get_matches_table()",
            )))
        }
    }

    /// Get similar table.
    pub fn get_similar_table(&self) -> Result<i32, AlignResultError> {
        if self.is_stats_table() {
            unsafe { Ok(*parasail_result_get_similar_table(self.inner)) }
        } else {
            Err(AlignResultError::NoStatsTable(String::from(
                "get_similar_table()",
            )))
        }
    }

    /// Get length table.
    pub fn get_length_table(&self) -> Result<i32, AlignResultError> {
        if self.is_stats_table() {
            unsafe { Ok(*parasail_result_get_length_table(self.inner)) }
        } else {
            Err(AlignResultError::NoStatsTable(String::from(
                "get_length_table()",
            )))
        }
    }

    /// Get score row.
    pub fn get_score_row(&self) -> Result<i32, AlignResultError> {
        if self.is_rowcol() || self.is_stats_rowcol() {
            unsafe { Ok(*parasail_result_get_score_row(self.inner)) }
        } else {
            Err(AlignResultError::NoRowCol(String::from("get_score_row()")))
        }
    }

    /// Get matches row.
    pub fn get_matches_row(&self) -> Result<i32, AlignResultError> {
        if self.is_stats_rowcol() {
            unsafe { Ok(*parasail_result_get_matches_row(self.inner)) }
        } else {
            Err(AlignResultError::NoRowCol(String::from(
                "get_matches_row()",
            )))
        }
    }

    /// Get similar row.
    pub fn get_similar_row(&self) -> Result<i32, AlignResultError> {
        if self.is_stats_rowcol() {
            unsafe { Ok(*parasail_result_get_similar_row(self.inner)) }
        } else {
            Err(AlignResultError::NoRowCol(String::from(
                "get_similar_row()",
            )))
        }
    }

    /// Get length row.
    pub fn get_length_row(&self) -> Result<i32, AlignResultError> {
        if self.is_stats_rowcol() {
            unsafe { Ok(*parasail_result_get_length_row(self.inner)) }
        } else {
            Err(AlignResultError::NoRowCol(String::from("get_length_row()")))
        }
    }

    /// Get score column.
    pub fn get_score_col(&self) -> Result<i32, AlignResultError> {
        if self.is_rowcol() || self.is_stats_rowcol() {
            unsafe { Ok(*parasail_result_get_score_col(self.inner)) }
        } else {
            Err(AlignResultError::NoRowCol(String::from("get_score_col()")))
        }
    }

    /// Get matches column.
    pub fn get_matches_col(&self) -> Result<i32, AlignResultError> {
        if self.is_stats_rowcol() {
            unsafe { Ok(*parasail_result_get_matches_col(self.inner)) }
        } else {
            Err(AlignResultError::NoRowCol(String::from(
                "get_matches_col()",
            )))
        }
    }

    /// Get similar column.
    pub fn get_similar_col(&self) -> Result<i32, AlignResultError> {
        if self.is_stats_rowcol() {
            unsafe { Ok(*parasail_result_get_similar_col(self.inner)) }
        } else {
            Err(AlignResultError::NoRowCol(String::from(
                "get_similar_col()",
            )))
        }
    }

    /// Get length column
    pub fn get_length_col(&self) -> Result<i32, AlignResultError> {
        if self.is_stats_rowcol() {
            unsafe { Ok(*parasail_result_get_length_col(self.inner)) }
        } else {
            Err(AlignResultError::NoRowCol(String::from("get_length_col()")))
        }
    }

    /// Get trace table.
    pub fn get_trace_table(&self) -> Result<i32, AlignResultError> {
        if self.is_trace() {
            unsafe { Ok(*parasail_result_get_trace_table(self.inner)) }
        } else {
            Err(AlignResultError::NoTrace(String::from("get_trace_table()")))
        }
    }

    /// Get alignment strings and statistics
    pub fn print_traceback(&self, query: &[u8], reference: &[u8]) {
        if self.is_trace() {
            let query_len = query.len() as i32;
            let ref_len = reference.len() as i32;
            let query = CString::new(query).unwrap();
            let reference = CString::new(reference).unwrap();
            let query_str = CString::new("Query:").unwrap();
            let ref_str = CString::new("Target:").unwrap();
            let match_char = CString::new("|").unwrap();
            let mismatch_char = CString::new(" ").unwrap();
            let width = 80;
            let name_width = 7;
            let use_stats = 1;
            unsafe {
                parasail_traceback_generic(
                    query.as_ptr(),
                    query_len,
                    reference.as_ptr(),
                    ref_len,
                    query_str.as_ptr(),
                    ref_str.as_ptr(),
                    self.matrix,
                    self.inner,
                    *match_char.as_ptr(),
                    *mismatch_char.as_ptr(),
                    *mismatch_char.as_ptr(),
                    width,
                    name_width,
                    use_stats,
                );
            }
        } else {
            println!("Alignment string is not available without traceback enabled. Consider using the `use_trace` method on AlignerBuilder.");
        }
    }

    /// Get alignment strings.
    pub fn get_traceback_strings(
        &self,
        query: &[u8],
        reference: &[u8],
    ) -> Result<Traceback, AlignResultError> {
        if self.is_trace() {
            let query_len = query.len() as i32;
            let ref_len = reference.len() as i32;
            let query = CString::new(query).unwrap();
            let reference = CString::new(reference).unwrap();
            let match_char = CString::new("|").unwrap();
            let mismatch_char = CString::new(" ").unwrap();
            unsafe {
                let alignment = parasail_result_get_traceback(
                    self.inner,
                    query.as_ptr(),
                    query_len,
                    reference.as_ptr(),
                    ref_len,
                    self.matrix,
                    *match_char.as_ptr(),
                    *mismatch_char.as_ptr(),
                    *mismatch_char.as_ptr(),
                );

                let query_traceback = CString::from_raw((*alignment).query).into_string().unwrap();
                let comparison_traceback =
                    CString::from_raw((*alignment).comp).into_string().unwrap();
                let reference_traceback =
                    CString::from_raw((*alignment).ref_).into_string().unwrap();

                Ok(Traceback {
                    query: query_traceback.clone(),
                    comparison: comparison_traceback.clone(),
                    reference: reference_traceback.clone(),
                })
            }
        } else {
            return Err(AlignResultError::NoTrace(String::from(
                "get_traceback_strings()",
            )));
        }
    }

    /// Get CIGAR string.
    pub fn get_cigar(&self, query: &[u8], reference: &[u8]) -> Result<String, AlignResultError> {
        if self.is_trace() {
            let query_len = query.len() as i32;
            let query = CString::new(query)?;
            let ref_len = reference.len() as i32;
            let reference = CString::new(reference)?;

            let cigar: Result<String, IntoStringError>;
            unsafe {
                let cigar_encoded = CigarString {
                    inner: parasail_result_get_cigar(
                        self.inner,
                        query.as_ptr(),
                        query_len,
                        reference.as_ptr(),
                        ref_len,
                        self.matrix,
                    ),
                };

                cigar = CString::from_raw(parasail_cigar_decode(cigar_encoded.inner)).into_string();
            }

            Ok(cigar?)
        } else {
            Err(AlignResultError::NoTrace(String::from("get_cigar()")))
        }
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
impl Drop for AlignResult {
    fn drop(&mut self) {
        unsafe {
            parasail_result_free(self.inner);
        }
    }
}

/// SSW alignment result.
pub struct SSWResult {
    inner: *mut parasail_result_ssw_t,
}

impl SSWResult {
    /// Get primary alignment score.
    pub fn score(&self) -> u16 {
        unsafe { (*self.inner).score1 }
    }

    /// Get beginning location of alignment on the reference sequence.
    pub fn ref_start(&self) -> i32 {
        unsafe { (*self.inner).ref_begin1 }
    }

    /// Get ending location of alignment on the reference sequence.
    pub fn ref_end(&self) -> i32 {
        unsafe { (*self.inner).ref_end1 }
    }

    /// Get beginning location of alignment on the query sequence.
    pub fn query_start(&self) -> i32 {
        unsafe { (*self.inner).read_begin1 }
    }

    /// Get ending location of alignment on the query sequence.
    pub fn query_end(&self) -> i32 {
        unsafe { (*self.inner).read_end1 }
    }

    //pub fn cigar(&self) -> *mut u32 {
    //    unsafe { (*self.inner).cigar }
    //}

    //pub fn cigar_len(&self) -> i32 {
    //    unsafe { (*self.inner).cigarLen }
    //}
}

#[doc(hidden)]
impl Drop for SSWResult {
    fn drop(&mut self) {
        unsafe { parasail_result_ssw_free(self.inner) }
    }
}
