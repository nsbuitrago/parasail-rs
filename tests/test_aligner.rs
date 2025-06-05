extern crate parasail_rs;

use std::thread;

use parasail_rs::{Aligner, Matrix, Profile};

const DNA_QUERY: &[u8; 5] = b"ACGGT";
const DNA_REFERENCE: &[u8; 5] = b"ACGGT";

const PROTEIN_QUERY: &[u8; 7] = b"MAGICAI"; // I -> L should be tolerated (similar) with Blosum62
const PROTEIN_REFERENCE: &[u8; 7] = b"MAGICAL";

#[test]
fn build_default_aligner() {
    assert!(Aligner::new().build().is_ok()); // default is global aligner
    assert!(Aligner::new().local().build().is_ok());
    assert!(Aligner::new().semi_global().build().is_ok());
}

#[test]
fn default_alignment() -> Result<(), Box<dyn std::error::Error>> {
    let global_aligner = Aligner::new().build()?;
    let local_aligner = Aligner::new().local().build()?;
    let semi_global_aligner = Aligner::new().semi_global().build()?;

    let global_result = global_aligner.align(DNA_QUERY, DNA_REFERENCE)?;
    let local_result = local_aligner.align(DNA_QUERY, DNA_REFERENCE)?;
    let semi_global_result = semi_global_aligner.align(DNA_QUERY, DNA_REFERENCE)?;

    assert!(global_result.is_global());
    assert!(!global_result.is_local());
    assert!(!global_result.is_semi_global());

    assert!(local_result.is_local());
    assert!(!local_result.is_global());
    assert!(!local_result.is_semi_global());

    assert!(semi_global_result.is_semi_global());
    assert!(!semi_global_result.is_local());
    assert!(!semi_global_result.is_global());

    Ok(())
}

#[test]
fn get_alignment_score() -> Result<(), Box<dyn std::error::Error>> {
    let aligner = Aligner::new().build()?;
    let alignment = aligner.align(DNA_QUERY, DNA_REFERENCE)?;

    assert_eq!(alignment.get_score(), DNA_QUERY.len() as i32);

    Ok(())
}

#[test]
fn get_end_query() -> Result<(), Box<dyn std::error::Error>> {
    let aligner = Aligner::new().build()?;
    let alignment = aligner.align(DNA_QUERY, DNA_REFERENCE)?;

    assert_eq!(alignment.get_end_query(), (DNA_QUERY.len() - 1) as i32);

    Ok(())
}

#[test]
fn get_end_ref() -> Result<(), Box<dyn std::error::Error>> {
    let aligner = Aligner::new().build()?;
    let alignment = aligner.align(DNA_QUERY, DNA_REFERENCE)?;

    assert_eq!(alignment.get_end_ref(), (DNA_REFERENCE.len() - 1) as i32);

    Ok(())
}

#[test]
fn get_matches() -> Result<(), Box<dyn std::error::Error>> {
    let aligner = Aligner::new().build()?;
    let default_alignment = aligner.align(DNA_QUERY, DNA_REFERENCE)?;
    let failed_n_matches = default_alignment.get_matches();

    let aligner_w_stats = Aligner::new().return_stats().build()?;
    let alignment = aligner_w_stats.align(DNA_QUERY, DNA_REFERENCE)?;
    let n_matches = alignment.get_matches();

    assert!(failed_n_matches.is_err());
    assert!(n_matches.is_ok());
    assert_eq!(n_matches.unwrap(), DNA_QUERY.len() as i32);

    Ok(())
}

#[test]
fn get_similar() -> Result<(), Box<dyn std::error::Error>> {
    let dna_aligner = Aligner::new().return_stats().build()?;
    let dna_alignment = dna_aligner.align(DNA_QUERY, DNA_REFERENCE)?;

    let invalid_test_case = dna_alignment.get_similar();
    println!("{}", invalid_test_case.as_ref().unwrap());
    assert!(invalid_test_case.is_ok());

    let blosum_matrix = Matrix::from("blosum62")?;
    let prot_aligner = Aligner::new()
        .matrix(blosum_matrix)
        .return_stats()
        .build()?;
    let prot_alignment = prot_aligner.align(PROTEIN_QUERY, PROTEIN_REFERENCE)?;
    let n_similar_residues = prot_alignment.get_similar();

    assert!(n_similar_residues.is_ok());
    assert_eq!(n_similar_residues.unwrap(), PROTEIN_QUERY.len() as i32);

    Ok(())
}

#[test]
fn get_length() -> Result<(), Box<dyn std::error::Error>> {
    let dna_aligner = Aligner::new().return_stats().build()?;
    let dna_alignment = dna_aligner.align(DNA_QUERY, DNA_REFERENCE)?;
    let length = dna_alignment.get_length();

    assert!(length.is_ok());
    assert_eq!(length.unwrap(), DNA_QUERY.len() as i32);
    Ok(())
}

#[test]
fn get_score_table() -> Result<(), Box<dyn std::error::Error>> {
    let aligner = Aligner::new().return_stats().build()?;
    let alignment = aligner.align(DNA_QUERY, DNA_REFERENCE)?;

    let score_table = alignment.get_score_table();

    assert!(score_table.is_ok());
    println!("{}", score_table.unwrap());

    Ok(())
}

#[test]
fn print_traceback() -> Result<(), Box<dyn std::error::Error>> {
    let aligner = Aligner::new().enable_traceback().build()?;
    let alignment = aligner.align(DNA_QUERY, DNA_REFERENCE)?;

    // alignment.print_traceback()?;
    //

    Ok(())
}

#[test]
fn profile_alignment() -> Result<(), Box<dyn std::error::Error>> {
    let query = b"ACTTAGCTT";
    let reference = b"ACTTAGCTT";
    let matrix = Matrix::default();
    let query_profile = Profile::new(query, &matrix).build()?;
    let query_profile_2 = query_profile.clone();
    let query_profile_3 = query_profile.clone();

    // alignment with profile requires scan or striped methods
    let global_aligner = Aligner::new().profile(query_profile).scan().build()?;

    let semi_global_aligner = Aligner::new()
        .semi_global()
        .profile(query_profile_2)
        .scan()
        .build()?;

    let local_aligner = Aligner::new()
        .local()
        .profile(query_profile_3)
        .scan()
        .build()?;

    let global_result = global_aligner.align_with_profile(reference)?;
    let semi_global_result = semi_global_aligner.align_with_profile(reference)?;
    let local_result = local_aligner.align_with_profile(reference)?;

    let target_score = query.len() as i32;
    assert_eq!(global_result.get_score(), target_score);
    assert_eq!(semi_global_result.get_score(), target_score);
    assert_eq!(local_result.get_score(), target_score);

    Ok(())
}

#[test]
pub fn multithreaded_global_alignment() -> Result<(), Box<dyn std::error::Error>> {
    let query_seq = b"ACGTGCA";
    let ref_seqs = vec![b"ACGTGCA", b"ACGTCGA"];

    let matrix = Matrix::default();
    let global_profile = Profile::new(query_seq, &matrix).use_stats().build()?;
    let semi_global_profile = Profile::new(query_seq, &matrix).use_stats().build()?;
    let local_profile = Profile::new(query_seq, &matrix).use_stats().build()?;

    let global_aligner = Aligner::new().profile(global_profile).striped().build()?;

    let semi_global_aligner = Aligner::new()
        .semi_global()
        .profile(semi_global_profile)
        .striped()
        .build()?;

    let local_aligner = Aligner::new()
        .local()
        .profile(local_profile)
        .striped()
        .build()?;

    thread::spawn(move || {
        for reference in ref_seqs {
            let global_result = global_aligner.align_with_profile(reference);
            let semi_global_result = semi_global_aligner.align_with_profile(reference);
            let local_result = local_aligner.align_with_profile(reference);

            assert!(global_result.is_ok());
            assert!(semi_global_result.is_ok());
            assert!(local_result.is_ok());
        }
    })
    .join()
    .unwrap();

    Ok(())
}
