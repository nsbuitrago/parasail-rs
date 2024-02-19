//! # Introduction
//!
//! This crate provides safe Rust bindings to the [Parasail](https://github.com/jeffdaily/parasail), a SIMD C library for //! pairwise sequence alignments.
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
//! TODO: update examples
//! ```rust,no_run
//! use parasail_rs::{ScoringMatrix, Aligner};
//! let matrix = ScoringMatrix::create("ACGT", 1, -1);
//! let query = b"ACGT";
//! let target = b"ACGTAACGTACA";
//!
//! let aligner = Aligner::new(matrix, 5, 2);
//! let result = aligner.global(query, target, "striped");
//! ```

use libparasail_sys::{
    parasail_matrix_create,
    parasail_matrix_lookup,
    parasail_matrix_from_file,
    parasail_matrix_pssm_create,
    parasail_matrix_convert_square_to_pssm,
    parasail_matrix_free,
    parasail_matrix_t,
    parasail_result_t,
    parasail_nw_striped_sat,
    parasail_nw_scan_sat,
    parasail_result_get_score,
    parasail_result_get_end_query,
    parasail_result_get_end_ref,
    parasail_result_get_matches,
    parasail_result_free,
    parasail_result_is_diag,
    parasail_result_is_saturated,
    parasail_result_is_sg,
    parasail_result_is_banded,
    parasail_result_is_striped,
    parasail_result_is_scan,
    parasail_result_is_nw,
    parasail_result_is_sw,
    parasail_result_get_length,
    parasail_result_get_similar,
};

#[derive(Debug)]
pub struct ScoringMatrix {
    pub scores: *mut parasail_matrix_t,
}

impl ScoringMatrix {
    pub fn create(alphabet: &str, match_score: i32, mismatch_score: i32) -> ScoringMatrix {
        unsafe {
            let scores = parasail_matrix_create(alphabet.as_ptr() as *const i8, match_score, mismatch_score);
            ScoringMatrix { scores }
        }
    }

    pub fn from(matrix_name: &str) -> ScoringMatrix {
        unsafe {
            let scores = parasail_matrix_lookup(matrix_name.as_ptr() as *const i8);
            ScoringMatrix { scores: scores as *mut libparasail_sys::parasail_matrix_t}
        }
    }

    pub fn from_file(file: &str) -> ScoringMatrix {
        unsafe {
            let scores = parasail_matrix_from_file(file.as_ptr() as *const i8);
            ScoringMatrix { scores }
        }
    }

    pub fn create_pssm(alphabet: &str, values: Vec<i32>, rows: i32) -> ScoringMatrix {
        unsafe {
            let scores = parasail_matrix_pssm_create(alphabet.as_ptr() as *const i8, values.as_ptr(), rows);
            ScoringMatrix { scores }
        }
    }

    pub fn convert_square_to_pssm(&self, pssm_query: &str) -> ScoringMatrix {
        unsafe {
            let scores = parasail_matrix_convert_square_to_pssm(self.scores, pssm_query.as_ptr() as *const i8, pssm_query.len() as i32);
            ScoringMatrix { scores }
        }
    }

}

impl Drop for ScoringMatrix {
    fn drop(&mut self) {
        unsafe {
            parasail_matrix_free(self.scores);
        }
    }
}

pub struct Aligner {
    pub matrix: ScoringMatrix,
    pub gap_open: i32,
    pub gap_extend: i32,
}


impl Aligner {
    pub fn new(matrix: ScoringMatrix, gap_open: i32, gap_extend: i32) -> Aligner {
        Aligner { matrix, gap_open, gap_extend }
    }

    pub fn global(&self, query: &[u8], target: &[u8], method: &str) -> AlignmentResult {
        let query_len = query.len() as i32;
        let target_len = target.len() as i32;
        unsafe {
            match method {
                "striped" => {
                    let result = parasail_nw_striped_sat(query.as_ptr() as *const i8, query_len, target.as_ptr() as *const i8, target_len, self.gap_open, self.gap_extend, self.matrix.scores);
                    AlignmentResult { result }
                },
                "scan" => {
                    let result = parasail_nw_scan_sat(query.as_ptr() as *const i8, query_len, target.as_ptr() as *const i8, target_len, self.gap_open, self.gap_extend, self.matrix.scores);
                    AlignmentResult { result }
                },
                _ => {
                    panic!("Invalid alignment method");
                }
                
            }
        }
    }
}

#[derive(Debug)]
pub struct AlignmentResult {
    pub result: *mut parasail_result_t,
}

impl AlignmentResult {
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

    pub fn is_global(&self) -> bool {
        unsafe {
            parasail_result_is_nw(self.result) == 1
        }
    }

    pub fn is_semi_global(&self) -> bool {
        unsafe {
            parasail_result_is_sg(self.result) == 1
        }
    }

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

    pub fn is_scan(&self) -> bool {
        unsafe {
            parasail_result_is_scan(self.result) == 1
        }
    }

    pub fn is_striped(&self) -> bool {
        unsafe {
            parasail_result_is_striped(self.result) == 1
        }
    }

    pub fn is_diag(&self) -> bool {
        unsafe {
            parasail_result_is_diag(self.result) == 1
        }
    }
}

impl Drop for AlignmentResult {
    fn drop(&mut self) {
        unsafe {
            parasail_result_free(self.result);
        }
    }
}

