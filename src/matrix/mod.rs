mod error;

pub use error::*; // flatten
use libparasail_sys::{
    parasail_matrix_copy, parasail_matrix_create, parasail_matrix_free, parasail_matrix_from_file,
    parasail_matrix_lookup, parasail_matrix_pssm_create, parasail_matrix_set_value,
    parasail_matrix_t,
};
use std::{
    ffi::{c_int, CString},
    fmt::Display,
    ops::Deref,
    path::Path,
};

/// Substitution matrix for alignments. Several built-in matrices (i.e., blosum62) are available and
/// can be constructed with the `from` method. Alternatively, build a custom matrix with the `create`
/// method.
///
/// Examples:
/// ```rust,no_run
/// use parasail_rs::Matrix;
///
///  // use a built-in matrix
///  let blosum62_matrix = Matrix::from("blosum62")?;
///  println!(":blosum62_matrix?");
///
///  // create a custom matrix
///  let alphabet = b"ACGT";
///  let match_score = 3;
///  let mismatch_score = -2;
///  let custom_matrix = Matrix::create(alphabet, match_score, mismatch_score)?;
///  println!(":custom_matrix?");
///
///  # Ok::<(), parasail_rs::Error>(())
/// ```
#[non_exhaustive]
#[derive(Debug)]
pub enum Matrix {
    Builtin(*const parasail_matrix_t),
    Custom(*mut parasail_matrix_t),
}

impl Matrix {
    /// Create a new custom substitution matrix from a given an alphabet, match, and mismatch score.
    /// This method will return an error if the alphabet has internal null bytes or
    /// `libparasail_sys::parasail_matrix_create` returns a null pointer.
    pub fn create(alphabet: &[u8], match_score: i32, mismatch_score: i32) -> Result<Self> {
        let alphabet_c = CString::new(alphabet)?;
        let matrix = unsafe {
            parasail_matrix_create(
                alphabet_c.as_ptr(),
                match_score as c_int,
                mismatch_score as c_int,
            )
        };
        if matrix.is_null() {
            Err(Error::MatrixCreationFailed {
                alphabet: String::from_utf8_lossy(alphabet).into(),
                match_score,
                mismatch_score,
            })
        } else {
            Ok(Matrix::Custom(matrix))
        }
    }

    /// Create a new position-specific scoring matrix from a given alphabet,
    /// ...
    pub fn create_pssm(alphabet: &[u8], values: Vec<i32>, n_rows: i32) -> Result<Self> {
        let pssm_alphabet = CString::new(alphabet)?;
        let pssm_matrix = unsafe {
            parasail_matrix_pssm_create(pssm_alphabet.as_ptr(), values.as_ptr(), n_rows as c_int)
        };

        if pssm_matrix.is_null() {
            Err(Error::PSSMCreationFailed {
                alphabet: String::from_utf8_lossy(alphabet).into(),
                values,
                n_rows,
            })
        } else {
            Ok(Matrix::Custom(pssm_matrix))
        }
    }

    /// Create a new built-in substitution matrix. This method will return an error if any of the
    /// below conditions are met:
    /// - a matrix with <name> does not exist
    /// - the matrix name has internal null bytes
    /// - `libparasail_sys::parasail_matrix_lookup` returns a null pointer.
    ///
    /// Below is a full list of available built-in matrices:
    /// - blosum30, blosum35, ..., blosum62, blosum65, blosum70, ..., blosum90
    /// - pam10, pam20, ..., pam500
    /// - dnafull
    /// - nuc44
    pub fn from(name: &str) -> Result<Self> {
        let name_c = CString::new(name)?;
        let matrix = unsafe { parasail_matrix_lookup(name_c.as_ptr()) };
        if matrix.is_null() {
            Err(Error::MatrixLookupFailed {
                matrix_name: name.to_string(),
            })
        } else {
            Ok(Matrix::Builtin(matrix))
        }
    }

    /// Create a square to position-specific scoring matrix from a file.
    pub fn from_file<P: AsRef<Path>>(filepath: P) -> Result<Self> {
        let filepath_c = CString::new(filepath.as_ref().display().to_string())?;
        let matrix = unsafe { parasail_matrix_from_file(filepath_c.as_ptr()) };

        if matrix.is_null() {
            Err(Error::FromFileFailed {
                filepath: filepath.as_ref().display().to_string(),
            })
        } else {
            Ok(Matrix::Custom(matrix))
        }
    }

    /// Set the value of a matrix at a given row and column index.
    pub fn set_value(&mut self, row: i32, col: i32, value: i32) -> Result<()> {
        if let Matrix::Builtin(_) = self {
            *self = self.clone()
        }

        let parasail_matrix = **self as *mut parasail_matrix_t;
        let size = unsafe { (*parasail_matrix).size - 2 };
        indices_are_valid(size, row, col)?;
        unsafe { parasail_matrix_set_value(parasail_matrix, row, col, value) }

        Ok(())
    }

    /// Get value of a matrix at a given row and column index.
    pub fn get_value(&self, row: i32, col: i32) -> Result<i32> {
        let parasail_matrix = **self as *mut parasail_matrix_t;
        let size = unsafe { (*parasail_matrix).size - 2 };
        indices_are_valid(size, row, col)?;

        let value = unsafe {
            let value_ptr = (*parasail_matrix)
                .matrix
                .offset((row * (size + 1) + col) as isize);
            *value_ptr
        };

        Ok(value)
    }
}

fn indices_are_valid(size: i32, row: i32, col: i32) -> Result<()> {
    if row < 0 || row > size || col < 0 || col > size {
        return Err(Error::InvalidIndex { row, col });
    }

    Ok(())
}

impl Default for Matrix {
    /// Creates an identity matrix for the DNA alphabet.
    fn default() -> Self {
        Matrix::create(b"ACGT", 1, -1).unwrap() // is there a better way to do this?
    }
}

#[doc(hidden)]
impl Display for Matrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        unsafe {
            let parasail_matrix = **self;
            let length = (*parasail_matrix).length as usize;
            let size = (*parasail_matrix).size as usize;
            let matrix = std::slice::from_raw_parts((*parasail_matrix).matrix, length * size);
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
        match self {
            Matrix::Custom(m) => unsafe {
                std::mem::transmute::<&*mut parasail_matrix_t, &*const parasail_matrix_t>(m)
            },
            Matrix::Builtin(m) => m,
        }
    }
}

#[doc(hidden)]
impl Clone for Matrix {
    fn clone(&self) -> Self {
        let parasail_matrix_copy = unsafe { parasail_matrix_copy(**self) };
        match *self {
            Matrix::Builtin(_) => Matrix::Builtin(parasail_matrix_copy as *const parasail_matrix_t),
            Matrix::Custom(_) => Matrix::Custom(parasail_matrix_copy),
        }
    }
}

#[doc(hidden)]
impl Drop for Matrix {
    // built-in matrices do not need to be freed as they are not allocated on the heap
    fn drop(&mut self) {
        if let Matrix::Custom(matrix) = *self {
            unsafe { parasail_matrix_free(matrix) }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[ignore]
    fn display_matrix() {
        let matrix = Matrix::default();
        println!("{matrix}");
    }
}
