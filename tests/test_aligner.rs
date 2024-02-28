use parasail_rs::{Matrix, Profile, Aligner};
use std::thread;

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

#[test]
pub fn test_multithread_alignment() {
    let query = b"ACGT";
    let refs = vec![b"ACGT", b"ACGT"];
    let matrix = Matrix::default();
    let profile = Profile::new(query, false, &matrix);

    let aligner = Aligner::new()
        .profile(profile)
        .build();
    
    thread::spawn(move || {
        for reference in refs {
            let result = &aligner.global_with_profile(reference);
            let score = result.get_score();
            assert_eq!(score, query.len() as i32);
        }
    }).join().unwrap();
}
