use parasail_rs::{Matrix, Profile, Aligner};
use std::io;

#[test]
pub fn test_matrix_construction() {
    let default_matrix = Matrix::default();
    //todo!()
}

#[test]
pub fn test_profile_construction() {
    let query = b"ATGGCACTATAA";
    let use_stats = false;
    let profile = Profile::new(query, use_stats, &Matrix::default());
    // todo!()
}

#[test]
pub fn test_aligner_construction() {
    let aligner = Aligner::new().build();
    //todo!()
}

#[test]
pub fn test_global_alignment() -> Result<(), io::Error> {

    let matrix = Matrix::default();
    let query = b"ACGT";
    let reference = b"ACGT";
    let with_stats = false;
    let profile = Profile::new(query, with_stats, &matrix);

    let aligner = Aligner::new().build();

    aligner.global(query, reference);

    Ok(())
}
