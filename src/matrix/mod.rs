mod error;

use crate::Result;
pub use error::Error;
use libparasail_sys::{
    parasail_matrix_convert_square_to_pssm, parasail_matrix_copy, parasail_matrix_create,
    parasail_matrix_free, parasail_matrix_from_file, parasail_matrix_lookup,
    parasail_matrix_pssm_create, parasail_matrix_set_value, parasail_matrix_t,
};
use std::ffi::CString;
use std::fmt::Display;
use std::ops::Deref;
use std::path::Path;
use std::slice;

/// Substitution matrix for sequence alignment.
/// Matrices can be created from:
/// - an alphabet and match/mismatch scores
/// - a pre-defined matrix (such as blosum62)
/// - a file containing a substitution matrix (see `from_file` for details)
/// - a PSSM (position-specific scoring matrix)
#[derive(Debug)]
pub struct Matrix {
    pub(crate) inner: *const parasail_matrix_t,
    pub(crate) builtin: bool,
}

impl Matrix {
    /// Create a new scoring matrix from an alphabet and match/mismatch scores.
    /// Note that match score should be a positive integer, while mismatch score
    /// should be a negative integer.
    pub fn create(alphabet: &[u8], match_score: i32, mismatch_score: i32) -> Result<Self> {
        assert!(match_score >= 0 && mismatch_score <= 0, "Match score should be a positive integer and mismatch score should be a negative integer.");
        assert!(!alphabet.is_empty(), "Alphabet should not be empty.");
        unsafe {
            let alphabet = &CString::new(alphabet).map_err(Error::CreateErr)?;
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
    pub fn from(matrix_name: &str) -> Result<Self> {
        assert!(!matrix_name.is_empty(), "Matrix name should not be empty.");
        let matrix: *const parasail_matrix_t;
        unsafe {
            let matrix_name = CString::new(matrix_name).map_err(Error::CreateErr)?;
            matrix = parasail_matrix_lookup(matrix_name.as_ptr());
        }

        if matrix.is_null() {
            return Err(Error::LookupErr(matrix_name.to_string()).into());
        }

        Ok(Self {
            inner: matrix,
            builtin: true,
        })
    }

    /// Create a new scoring matrix from a file.
    /// Files should contain either square or position-specific scoring matrices.
    /// The examples below are taken directly from the [Parasail C library docs](https://github.com/jeffdaily/parasail?tab=readme-ov-file#substitution-matrices).
    ///
    /// A square matrix can be loaded from a file with the following contents:
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
    /// A position specific scoring matrix can also be loaded from a file.
    /// For example:
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
    pub fn from_file(file: &str) -> Result<Self> {
        let filepath = Path::new(file);
        if !filepath.exists() {
            return Err(Error::FileNotFound(filepath.to_str().unwrap_or("").to_string()).into());
        }

        let file = CString::new(file).map_err(Error::CreateErr)?;

        unsafe {
            let matrix = parasail_matrix_from_file(file.as_ptr());

            if matrix.is_null() {
                return Err(Error::FromFileErr(filepath.to_str().unwrap_or("").to_string()).into());
            }

            Ok(Self {
                inner: parasail_matrix_from_file(file.as_ptr()),
                builtin: false,
            })
        }
    }

    /// Create a new scoring matrix from a position-specific scoring matrix.
    pub fn create_pssm(alphabet: &str, values: Vec<i32>, rows: i32) -> Result<Self> {
        let alphabet = CString::new(alphabet).map_err(Error::CreateErr)?;

        unsafe {
            let matrix = parasail_matrix_pssm_create(alphabet.as_ptr(), values.as_ptr(), rows);

            if matrix.is_null() {
                return Err(Error::CreatePssmErr.into());
            }

            Ok(Self {
                inner: matrix,
                builtin: false,
            })
        }
    }

    /// Convert a square scoring matrix to a PSSM (position-specific scoring matrix).
    pub fn to_pssm(self, pssm_query: &[u8]) -> Result<Matrix> {
        assert!(
            !pssm_query.is_empty(),
            "PSSM query sequence should not be empty."
        );
        let pssm_query_string = CString::new(pssm_query).map_err(Error::CreateErr)?;

        unsafe {
            let matrix = parasail_matrix_copy(self.inner);
            if matrix.is_null() {
                return Err(Error::NullMatrix.into());
            }

            if (*self.inner).type_ != 0 {
                return Err(Error::NotSquare.into());
            }

            let converted_matrix = parasail_matrix_convert_square_to_pssm(
                matrix,
                pssm_query_string.as_ptr(),
                pssm_query.len() as i32,
            );

            if converted_matrix.is_null() {
                panic!("Error converting matrix to PSSM. Invalid query sequence.")
            }

            Ok(Matrix {
                inner: converted_matrix,
                builtin: self.builtin,
            })
        }
    }

    /// Set value at a given row and column index for a user defined substitution matrix.
    pub fn set_value(&mut self, row: i32, col: i32, value: i32) -> Result<()> {
        if self.builtin {
            return Err(Error::NotBuiltIn.into());
        }

        unsafe {
            let size = (*self.inner).size - 2;

            if size < 0 {
                return Err(Error::NullMatrix.into());
            }

            if row < 0 || row > size || col < 0 || col > size {
                return Err(Error::InvalidIndex(row, col, size).into());
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
impl Clone for Matrix {
    fn clone(&self) -> Self {
        let parasail_matrix_copy = unsafe { parasail_matrix_copy(**self) };
        if self.builtin {
            Self {
                inner: parasail_matrix_copy as *const parasail_matrix_t,
                builtin: false, // for consistency with C interface
            }
        } else {
            Self {
                inner: parasail_matrix_copy,
                builtin: false,
            }
        }
    }
}

#[doc(hidden)]
impl Drop for Matrix {
    fn drop(&mut self) {
        if self.builtin {
            return;
        }

        if !self.inner.is_null() {
            unsafe { parasail_matrix_free(self.inner as *mut parasail_matrix_t) }
        }
    }
}

#[doc(hidden)]
unsafe impl Send for Matrix {}
#[doc(hidden)]
unsafe impl Sync for Matrix {}
