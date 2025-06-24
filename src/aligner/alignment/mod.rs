mod error;

use libparasail_sys::{
    parasail_cigar_decode, parasail_cigar_free, parasail_cigar_t, parasail_matrix_t,
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
    parasail_traceback_generic,
};
use std::ffi::CString;
use std::slice;

use crate::Result;
pub use error::Error;

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
pub struct Alignment {
    pub(crate) inner: *mut parasail_result_t,
    pub(crate) matrix: *const parasail_matrix_t,
    pub(crate) query_len: i32,
    pub(crate) ref_len: i32,
}

impl Alignment {
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
    pub fn get_matches(&self) -> Result<i32> {
        if self.is_stats() {
            unsafe { Ok(parasail_result_get_matches(self.inner)) }
        } else {
            Err(Error::NoStats(String::from("get_matches()")))?
        }
    }

    pub fn get_similar(&self) -> i32 {
        unsafe { parasail_result_get_similar(self.inner) }
    }

    /// Get alignment length.
    pub fn get_length(&self) -> Result<i32> {
        if self.is_stats() {
            unsafe { Ok(parasail_result_get_length(self.inner)) }
        } else {
            Err(Error::NoStats(String::from("get_length()")))?
        }
    }

    /// Get score table.
    pub fn get_score_table(&self) -> Result<i32> {
        if self.is_table() || self.is_stats_table() {
            unsafe { Ok(*parasail_result_get_score_table(self.inner)) }
        } else {
            Err(Error::NoTable(String::from("get_score_table()")))?
        }
    }

    /// Get matches table.
    pub fn get_matches_table(&self) -> Result<i32> {
        if self.is_stats_table() {
            unsafe { Ok(*parasail_result_get_matches_table(self.inner)) }
        } else {
            Err(Error::NoStatsTable(String::from("get_matches_table()")))?
        }
    }

    /// Get similar table.
    pub fn get_similar_table(&self) -> Result<i32> {
        if self.is_stats_table() {
            unsafe { Ok(*parasail_result_get_similar_table(self.inner)) }
        } else {
            Err(Error::NoStatsTable(String::from("get_similar_table()")))?
        }
    }

    /// Get length table.
    pub fn get_length_table(&self) -> Result<i32> {
        if self.is_stats_table() {
            unsafe { Ok(*parasail_result_get_length_table(self.inner)) }
        } else {
            Err(Error::NoStatsTable(String::from("get_length_table()")))?
        }
    }

    /// Get score row.
    pub fn get_score_row(&self) -> Result<&[i32]> {
        if self.is_rowcol() || self.is_stats_rowcol() {
            unsafe {
                let rptr = parasail_result_get_score_row(self.inner);
                Ok(slice::from_raw_parts(rptr, self.ref_len as usize))
            }
        } else {
            Err(Error::NoRowCol(String::from("get_score_row()")))?
        }
    }

    /// Get matches row.
    pub fn get_matches_row(&self) -> Result<&[i32]> {
        if self.is_stats_rowcol() {
            unsafe {
                let rptr = parasail_result_get_matches_row(self.inner);
                Ok(slice::from_raw_parts(rptr, self.ref_len as usize))
            }
        } else {
            Err(Error::NoRowCol(String::from("get_matches_row()")))?
        }
    }

    /// Get similar row.
    pub fn get_similar_row(&self) -> Result<&[i32]> {
        if self.is_stats_rowcol() {
            unsafe {
                let rptr = parasail_result_get_similar_row(self.inner);
                Ok(slice::from_raw_parts(rptr, self.ref_len as usize))
            }
        } else {
            Err(Error::NoRowCol(String::from("get_similar_row()")))?
        }
    }

    /// Get length row.
    pub fn get_length_row(&self) -> Result<&[i32]> {
        if self.is_stats_rowcol() {
            unsafe {
                let rptr = parasail_result_get_length_row(self.inner);
                Ok(slice::from_raw_parts(rptr, self.ref_len as usize))
            }
        } else {
            Err(Error::NoRowCol(String::from("get_length_row()")))?
        }
    }

    /// Get score column.
    pub fn get_score_col(&self) -> Result<&[i32]> {
        if self.is_rowcol() || self.is_stats_rowcol() {
            unsafe {
                let cptr = parasail_result_get_score_col(self.inner);
                Ok(slice::from_raw_parts(cptr, self.query_len as usize))
            }
        } else {
            Err(Error::NoRowCol(String::from("get_score_col()")))?
        }
    }

    /// Get matches column.
    pub fn get_matches_col(&self) -> Result<&[i32]> {
        if self.is_stats_rowcol() {
            unsafe {
                let cptr = parasail_result_get_matches_col(self.inner);
                Ok(slice::from_raw_parts(cptr, self.query_len as usize))
            }
        } else {
            Err(Error::NoRowCol(String::from("get_matches_col()")))?
        }
    }

    /// Get similar column.
    pub fn get_similar_col(&self) -> Result<&[i32]> {
        if self.is_stats_rowcol() {
            unsafe {
                let cptr = parasail_result_get_similar_col(self.inner);
                Ok(slice::from_raw_parts(cptr, self.query_len as usize))
            }
        } else {
            Err(Error::NoRowCol(String::from("get_similar_col()")))?
        }
    }

    /// Get length column
    pub fn get_length_col(&self) -> Result<&[i32]> {
        if self.is_stats_rowcol() {
            unsafe {
                let cptr = parasail_result_get_length_col(self.inner);
                Ok(slice::from_raw_parts(cptr, self.query_len as usize))
            }
        } else {
            Err(Error::NoRowCol(String::from("get_length_col()")))?
        }
    }

    /// Get trace table.
    pub fn get_trace_table(&self) -> Result<i32> {
        if self.is_trace() {
            unsafe { Ok(*parasail_result_get_trace_table(self.inner)) }
        } else {
            Err(Error::NoTrace(String::from("get_trace_table()")))?
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
    pub fn get_traceback_strings(&self, query: &[u8], reference: &[u8]) -> Result<Traceback> {
        if self.is_trace() {
            let query_len = query.len() as i32;
            let ref_len = reference.len() as i32;
            let query = CString::new(query).map_err(Error::NewCStringErr)?;
            let reference = CString::new(reference).map_err(Error::NewCStringErr)?;
            let match_char = CString::new("|").map_err(Error::NewCStringErr)?;
            let mismatch_char = CString::new(" ").map_err(Error::NewCStringErr)?;
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

                let query_traceback = CString::from_raw((*alignment).query)
                    .into_string()
                    .map_err(Error::CigarToStringErr)?;
                let comparison_traceback = CString::from_raw((*alignment).comp)
                    .into_string()
                    .map_err(Error::CigarToStringErr)?;
                let reference_traceback = CString::from_raw((*alignment).ref_)
                    .into_string()
                    .map_err(Error::CigarToStringErr)?;

                Ok(Traceback {
                    query: query_traceback.clone(),
                    comparison: comparison_traceback.clone(),
                    reference: reference_traceback.clone(),
                })
            }
        } else {
            Err(Error::NoTrace(String::from("get_traceback_strings()")))?
        }
    }

    /// Get CIGAR string.
    pub fn get_cigar(&self, query: &[u8], reference: &[u8]) -> Result<String> {
        if self.is_trace() {
            let query_len = query.len() as i32;
            let query = CString::new(query).map_err(Error::NewCStringErr)?;
            let ref_len = reference.len() as i32;
            let reference = CString::new(reference).map_err(Error::NewCStringErr)?;

            let cigar: String;
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

                cigar = CString::from_raw(parasail_cigar_decode(cigar_encoded.inner))
                    .into_string()
                    .map_err(Error::CigarToStringErr)?;
            }

            Ok(cigar)
        } else {
            Err(Error::NoTrace(String::from("get_cigar()")).into())
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
impl Drop for Alignment {
    fn drop(&mut self) {
        unsafe {
            parasail_result_free(self.inner);
        }
    }
}

/// SSW alignment result.
pub struct SSWResult {
    pub(crate) inner: *mut parasail_result_ssw_t,
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

    pub fn cigar(&self) -> *mut u32 {
        unsafe { (*self.inner).cigar }
    }

    pub fn cigar_len(&self) -> i32 {
        unsafe { (*self.inner).cigarLen }
    }
}

#[doc(hidden)]
impl Drop for SSWResult {
    fn drop(&mut self) {
        unsafe { parasail_result_ssw_free(self.inner) }
    }
}
