use parasail_rs::{Matrix, Profile, Aligner};

#[test]
pub fn test_default_matrix_construction() {
    Matrix::default();
}

#[test]
pub fn test_profile_construction() {
    let query = b"ATGGCACTATAA";
    let use_stats = false;
    Profile::new(query, use_stats, &Matrix::default());
}

#[test]
pub fn test_aligner_construction() {
    Aligner::new().build();
}

#[test]
pub fn test_global_alignment() {
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().build();
    let result = aligner.global(query, reference);

    let score = result.get_score();
    assert_eq!(score, query.len() as i32);
}
