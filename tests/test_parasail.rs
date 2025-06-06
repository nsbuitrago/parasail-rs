use parasail_rs::{Aligner, Matrix, Profile};
use std::thread::{self};

#[test]
pub fn matrix_construction() -> Result<(), Box<dyn std::error::Error>> {
    // default matrix
    Matrix::default();

    // custom matrix
    let mut matrix = Matrix::create(b"ACGT", 3, -2)?;
    println!("Matrix:\n {}", matrix);
    matrix.set_value(2, 2, 100)?;
    println!("Matrix:\n {}", matrix);

    // from name
    let blosum62 = Matrix::from("blosum62")?;

    // convert to PSSM
    let _blosum62_pssm = blosum62.to_pssm(b"ACGT");

    // square matrix from file
    Matrix::from_file("./tests/square.txt")?;

    // pssm from file
    Matrix::from_file("./tests/pssm.txt")?;

    // PSSM
    let pssm_alphabet = "abcdef";
    let values = vec![1, 2, 3, 4, 5, 6, 7, 8];
    let rows = 2;
    Matrix::create_pssm(pssm_alphabet, values, rows)?;

    Ok(())
}

#[test]
pub fn profile_construction() -> Result<(), Box<dyn std::error::Error>> {
    let query = b"ATGGCACTATAA";
    Profile::new(query, false, &Matrix::default())?;

    // with stats
    Profile::new(query, true, &Matrix::default())?;

    Ok(())
}

#[test]
pub fn aligner_construction() {
    // default aligner
    Aligner::new().build();

    // custom aligner
    Aligner::new()
        .matrix(Matrix::default())
        .gap_open(10)
        .gap_extend(1)
        .profile(Profile::default())
        .allow_query_gaps(vec![String::from("prefix"), String::from("suffix")])
        .striped()
        .use_stats()
        .build();
}

#[test]
pub fn global_alignment() -> Result<(), Box<dyn std::error::Error>> {
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().striped().build();
    let result = aligner.align(Some(query), reference)?;

    let checks = query.len() as i32;

    assert_eq!(result.get_score(), checks);
    assert_eq!(result.get_end_query(), checks - 1);
    assert_eq!(result.get_end_ref(), checks - 1);
    assert!(result.is_global());
    assert!(!result.is_local());
    assert!(!result.is_semi_global());
    assert!(result.is_striped());

    Ok(())
}

#[test]
pub fn semi_global_alignment() -> Result<(), Box<dyn std::error::Error>> {
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().semi_global().striped().build();
    let result = aligner.align(Some(query), reference)?;

    let checks = query.len() as i32;

    assert_eq!(result.get_score(), checks);
    assert_eq!(result.get_end_query(), checks - 1);
    assert_eq!(result.get_end_ref(), checks - 1);
    assert!(!result.is_global());
    assert!(!result.is_local());
    assert!(result.is_semi_global());
    assert!(result.is_striped());

    Ok(())
}

#[test]
pub fn local_alignment() -> Result<(), Box<dyn std::error::Error>> {
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().local().striped().build();
    let result = aligner.align(Some(query), reference)?;

    let checks = query.len() as i32;

    assert_eq!(result.get_score(), checks);
    assert_eq!(result.get_end_query(), checks - 1);
    assert_eq!(result.get_end_ref(), checks - 1);
    assert!(!result.is_global());
    assert!(result.is_local());
    assert!(!result.is_semi_global());
    assert!(result.is_striped());

    Ok(())
}

#[test]
pub fn global_with_stats() -> Result<(), Box<dyn std::error::Error>> {
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().use_stats().striped().build();
    let result = aligner.align(Some(query), reference)?;

    let checks = query.len() as i32;
    let n_matches = result.get_matches()?;
    let align_len = result.get_length()?;

    assert_eq!(n_matches, checks);
    assert_eq!(align_len, checks);

    Ok(())
}

#[test]
pub fn semi_global_with_stats() -> Result<(), Box<dyn std::error::Error>> {
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().semi_global().use_stats().striped().build();
    let result = aligner.align(Some(query), reference)?;

    let checks = query.len() as i32;
    let n_matches = result.get_matches()?;
    let align_len = result.get_length()?;

    assert_eq!(n_matches, checks);
    assert_eq!(align_len, checks);

    Ok(())
}

#[test]
pub fn local_with_stats() -> Result<(), Box<dyn std::error::Error>> {
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().local().use_stats().striped().build();
    let result = aligner.align(Some(query), reference)?;

    let checks = query.len() as i32;
    let n_matches = result.get_matches()?;
    let align_len = result.get_length()?;

    assert_eq!(n_matches, checks);
    assert_eq!(align_len, checks);

    Ok(())
}

#[test]
pub fn global_8bit() -> Result<(), Box<dyn std::error::Error>> {
    let query = b"ACTGACTGACTG";
    let reference = b"ACTGTCTGACTG";
    let aligner = Aligner::new().striped().solution_width(8).build();
    let result = aligner.align(Some(query), reference)?;

    let checks = query.len() as i32;

    assert_eq!(result.get_score(), checks - 1);
    assert_eq!(result.get_end_query(), checks - 1);
    assert_eq!(result.get_end_ref(), checks - 1);
    assert!(result.is_global());
    assert!(!result.is_local());
    assert!(!result.is_semi_global());
    assert!(result.is_striped());

    Ok(())
}

#[test]
pub fn global_16bit() -> Result<(), Box<dyn std::error::Error>> {
    let query = b"ACTGACTGACTG";
    let reference = b"ACTGTCTGACTG";
    let aligner = Aligner::new().striped().solution_width(16).build();
    let result = aligner.align(Some(query), reference)?;

    let checks = query.len() as i32;

    assert_eq!(result.get_score(), checks - 1);
    assert_eq!(result.get_end_query(), checks - 1);
    assert_eq!(result.get_end_ref(), checks - 1);
    assert!(result.is_global());
    assert!(!result.is_local());
    assert!(!result.is_semi_global());
    assert!(result.is_striped());

    Ok(())
}

#[test]
pub fn global_32bit() -> Result<(), Box<dyn std::error::Error>> {
    let query = b"ACTGACTGACTG";
    let reference = b"ACTGTCTGACTG";
    let aligner = Aligner::new().striped().solution_width(32).build();
    let result = aligner.align(Some(query), reference)?;

    let checks = query.len() as i32;

    assert_eq!(result.get_score(), checks - 1);
    assert_eq!(result.get_end_query(), checks - 1);
    assert_eq!(result.get_end_ref(), checks - 1);
    assert!(result.is_global());
    assert!(!result.is_local());
    assert!(!result.is_semi_global());
    assert!(result.is_striped());

    Ok(())
}

#[test]
pub fn global_64bit() -> Result<(), Box<dyn std::error::Error>> {
    let query = b"ACTGACTGACTG";
    let reference = b"ACTGTCTGACTG";
    let aligner = Aligner::new().striped().solution_width(64).build();
    let result = aligner.align(Some(query), reference)?;

    let checks = query.len() as i32;

    assert_eq!(result.get_score(), checks - 1);
    assert_eq!(result.get_end_query(), checks - 1);
    assert_eq!(result.get_end_ref(), checks - 1);
    assert!(result.is_global());
    assert!(!result.is_local());
    assert!(!result.is_semi_global());
    assert!(result.is_striped());

    Ok(())
}

#[test]
pub fn score_table() -> Result<(), Box<dyn std::error::Error>> {
    // one-off alignment without stats
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().use_table().striped().build();
    let default_score = 1;

    let result = aligner.align(Some(query), reference)?;
    assert!(result.is_table());
    assert!(!result.is_stats());
    assert!(!result.is_stats_table());
    assert_eq!(result.get_score_table()?, default_score);

    // one-off alignment with stats
    let aligner = Aligner::new().use_stats().use_table().striped().build();

    let result = aligner.align(Some(query), reference)?;
    assert!(result.is_stats());
    assert!(result.is_stats_table());
    assert!(result.is_table());
    assert_eq!(result.get_score_table()?, default_score);

    // alignment with profile, without stats
    let custom_score = 3;
    let matrix = Matrix::create(b"ACGT", custom_score, -2)?;
    let profile = Profile::new(query, false, &matrix)?;

    let aligner_w_profile = Aligner::new()
        .profile(profile)
        .use_table()
        .striped()
        .build();

    let result_w_profile = aligner_w_profile.align(None, reference)?;

    assert!(result_w_profile.is_table());
    assert!(!result_w_profile.is_stats());
    assert!(!result_w_profile.is_stats_table());
    assert_eq!(result_w_profile.get_score_table()?, custom_score);

    // // alignment with profile, with stats
    let profile = Profile::new(query, true, &matrix)?;
    let aligner_w_profile = Aligner::new()
        .profile(profile)
        .use_stats()
        .use_table()
        .striped()
        .build();

    let result_w_profile = aligner_w_profile.align(None, reference)?;

    assert!(result_w_profile.is_stats());
    assert!(result_w_profile.is_stats_table());
    assert!(result_w_profile.is_table());
    assert_eq!(result_w_profile.get_score_table()?, custom_score);

    Ok(())
}

#[test]
pub fn matches_table() -> Result<(), Box<dyn std::error::Error>> {
    // one-off alignment without stats
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().use_table().use_stats().striped().build();
    let default_score = 1;

    let result = aligner.align(Some(query), reference)?;
    assert!(result.is_table());
    assert!(result.is_stats());
    assert!(result.is_stats_table());
    assert_eq!(result.get_matches_table()?, default_score);

    Ok(())
}

#[test]
pub fn similar_table() -> Result<(), Box<dyn std::error::Error>> {
    // one-off alignment without stats
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().use_table().use_stats().striped().build();

    let result = aligner.align(Some(query), reference)?;
    assert!(result.is_table());
    assert!(result.is_stats());
    assert!(result.is_stats_table());
    println!("similar table: {:?}", result.get_similar_table()?);

    Ok(())
}

#[test]
pub fn length_table() -> Result<(), Box<dyn std::error::Error>> {
    // one-off alignment without stats
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().use_table().use_stats().striped().build();

    let result = aligner.align(Some(query), reference)?;
    assert!(result.is_table());
    assert!(result.is_stats());
    assert!(result.is_stats_table());
    println!("Length table: {:?}", result.get_length_table()?);

    Ok(())
}

#[test]
pub fn score_row() -> Result<(), Box<dyn std::error::Error>> {
    // one-off alignment without stats
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new()
        .use_last_rowcol()
        .use_stats()
        .striped()
        .build();

    let result = aligner.align(Some(query), reference)?;
    assert!(result.is_stats_rowcol());
    assert!(result.is_stats());
    assert!(!result.is_stats_table());
    println!("Score row: {:?}", result.get_score_row()?);

    Ok(())
}

#[test]
pub fn matches_row() -> Result<(), Box<dyn std::error::Error>> {
    // one-off alignment without stats
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new()
        .use_last_rowcol()
        .use_stats()
        .striped()
        .build();

    let result = aligner.align(Some(query), reference)?;
    assert!(result.is_stats_rowcol());
    assert!(result.is_stats());
    assert!(!result.is_stats_table());
    println!("Matches row: {:?}", result.get_matches_row()?);

    Ok(())
}

#[test]
pub fn similar_row() -> Result<(), Box<dyn std::error::Error>> {
    // one-off alignment without stats
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new()
        .use_last_rowcol()
        .use_stats()
        .striped()
        .build();

    let result = aligner.align(Some(query), reference)?;
    assert!(result.is_stats_rowcol());
    assert!(result.is_stats());
    assert!(!result.is_stats_table());
    println!("Similar row: {:?}", result.get_similar_row()?);

    Ok(())
}

#[test]
pub fn length_row() -> Result<(), Box<dyn std::error::Error>> {
    // one-off alignment without stats
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new()
        .use_last_rowcol()
        .use_stats()
        .striped()
        .build();

    let result = aligner.align(Some(query), reference)?;
    assert!(result.is_stats_rowcol());
    assert!(result.is_stats());
    assert!(!result.is_stats_table());
    println!("Length row: {:?}", result.get_length_row());

    Ok(())
}

#[test]
pub fn score_col() -> Result<(), Box<dyn std::error::Error>> {
    // one-off alignment without stats
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new()
        .use_last_rowcol()
        .use_stats()
        .striped()
        .build();

    let result = aligner.align(Some(query), reference)?;
    assert!(result.is_stats_rowcol());
    assert!(result.is_stats());
    assert!(!result.is_stats_table());
    println!("Score col: {:?}", result.get_score_col()?);

    Ok(())
}

#[test]
pub fn match_col() -> Result<(), Box<dyn std::error::Error>> {
    // one-off alignment without stats
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new()
        .use_last_rowcol()
        .use_stats()
        .striped()
        .build();

    let result = aligner.align(Some(query), reference)?;
    assert!(result.is_stats_rowcol());
    assert!(result.is_stats());
    assert!(!result.is_stats_table());
    println!("Matches col: {:?}", result.get_matches_col()?);

    Ok(())
}

#[test]
pub fn similar_col() -> Result<(), Box<dyn std::error::Error>> {
    // one-off alignment without stats
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new()
        .use_last_rowcol()
        .use_stats()
        .striped()
        .build();

    let result = aligner.align(Some(query), reference)?;
    assert!(result.is_stats_rowcol());
    assert!(result.is_stats());
    assert!(!result.is_stats_table());
    println!("Similar col: {:?}", result.get_similar_col()?);

    Ok(())
}

#[test]
pub fn length_col() -> Result<(), Box<dyn std::error::Error>> {
    // one-off alignment without stats
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new()
        .use_last_rowcol()
        .use_stats()
        .striped()
        .build();

    let result = aligner.align(Some(query), reference)?;
    assert!(result.is_stats_rowcol());
    assert!(result.is_stats());
    assert!(!result.is_stats_table());
    println!("Length col: {:?}", result.get_length_col()?);

    Ok(())
}

#[test]
pub fn trace_table() -> Result<(), Box<dyn std::error::Error>> {
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().use_trace().striped().build();
    let result = aligner.align(Some(query), reference)?;
    assert!(result.is_trace());
    println!("Trace table: {:?}", result.get_trace_table()?);

    Ok(())
}

#[test]
pub fn get_traceback_strings() -> Result<(), Box<dyn std::error::Error>> {
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().use_trace().striped().build();
    let result = aligner.align(Some(query), reference)?;
    let traceback = result.get_traceback_strings(query, reference)?;
    println!("Query:     {}", traceback.query);
    println!("           {}", traceback.comparison);
    println!("Reference: {}", traceback.reference);

    Ok(())
}

#[test]
pub fn print_traceback() -> Result<(), Box<dyn std::error::Error>> {
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().use_trace().striped().build();
    let result = aligner.align(Some(query), reference)?;
    result.print_traceback(query, reference);

    Ok(())
}

#[test]
pub fn get_cigar() -> Result<(), Box<dyn std::error::Error>> {
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().use_trace().striped().build();
    let result = aligner.align(Some(query), reference)?;
    let cigar_string = result.get_cigar(query, reference)?;

    println!("CIGAR: {}", cigar_string);

    Ok(())
}

#[test]
pub fn global_with_profile() -> Result<(), Box<dyn std::error::Error>> {
    let query = b"ACGT";
    let reference = b"ACGT";
    let matrix = Matrix::default();
    let profile = Profile::new(query, true, &matrix)?;

    let aligner = Aligner::new()
        .profile(profile)
        .use_stats()
        .striped()
        .build();

    let result = aligner.align(None, reference)?;
    assert!(result.is_global());
    assert!(result.is_striped());
    assert!(result.is_stats());
    assert!(!result.is_local());
    assert!(!result.is_semi_global());

    Ok(())
}

#[test]
pub fn semi_global_with_profile() -> Result<(), Box<dyn std::error::Error>> {
    let query = b"ACGT";
    let reference = b"ACGT";
    let matrix = Matrix::default();
    let profile = Profile::new(query, true, &matrix)?;

    let aligner = Aligner::new()
        .profile(profile)
        .use_stats()
        .striped()
        .semi_global()
        .build();

    let result = aligner.align(None, reference)?;
    assert!(result.is_semi_global());
    assert!(result.is_striped());
    assert!(result.is_stats());
    assert!(!result.is_local());
    assert!(!result.is_global());

    Ok(())
}

#[test]
pub fn local_with_profile() -> Result<(), Box<dyn std::error::Error>> {
    let query = b"ACGT";
    let reference = b"ACGT";
    let matrix = Matrix::default();
    let profile = Profile::new(query, true, &matrix)?;

    let aligner = Aligner::new()
        .profile(profile)
        .use_stats()
        .striped()
        .local()
        .build();

    let result = aligner.align(None, reference)?;
    assert!(result.is_local());
    assert!(result.is_striped());
    assert!(result.is_stats());
    assert!(!result.is_global());
    assert!(!result.is_semi_global());

    Ok(())
}

#[test]
pub fn multithread_global_alignment() -> Result<(), Box<dyn std::error::Error>> {
    let query = b"ACGT";
    let refs = vec![b"ACGT", b"ACGT"];
    let matrix = Matrix::default();
    let profile = Profile::new(query, true, &matrix)?;

    let aligner = Aligner::new()
        .profile(profile)
        .use_stats()
        .scan()
        .build();

    let handles: Vec<_> = refs
        .iter()
        .map(|reference| {
            let aligner = aligner.clone();
            let reference = reference.to_vec();
            thread::spawn(move || match &aligner.align(None, &reference) {
                Ok(result) => {
                    let score = result.get_score();
                    assert_eq!(score, query.len() as i32);
                }
                Err(e) => {
                    println!("Alignment Error: {}", e);
                }
            })
        })
        .collect();

    for handle in handles {
        handle.join().unwrap();
    }
    Ok(())
}

#[test]
pub fn test_banded_nw() -> Result<(), Box<dyn std::error::Error>> {
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().bandwidth(2).build();
    let result = aligner.banded_nw(query, reference)?;
    let expected_score = query.len() as i32;

    assert_eq!(result.get_score(), expected_score);

    Ok(())
}

#[test]
pub fn test_ssw_alignment() -> Result<(), Box<dyn std::error::Error>> {
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().build();
    let result = aligner.ssw(Some(query), reference)?;

    let checks = query.len() as u16;
    let end = checks as i32 - 1;
    let start: i32 = 0;

    assert_eq!(result.score(), checks);
    assert_eq!(result.query_end(), end);
    assert_eq!(result.ref_end(), end);
    assert_eq!(result.query_start(), start);
    assert_eq!(result.ref_start(), start);

    Ok(())
}

#[test]
pub fn test_ssw_init() -> Result<(), Box<dyn std::error::Error>> {
    let query = b"ACGT";
    let matrix = Matrix::default();
    let _ = Profile::new_ssw(query, &matrix, 2);

    Ok(())
}

// #[test]
// pub fn test_ssw_profile_alignment() -> Result<(), Box<dyn std::error::Error>> {
//     let query = b"ACGTTGA";
//     let reference = b"ACGTTGA";
//     let matrix = Matrix::default();
//     let profile = Profile::new_ssw(query, &matrix, 1)?;
//
//     let aligner = Aligner::new().profile(profile).build();
//     let result = aligner.ssw(None, reference)?;
//
//     let checks = query.len() as u16;
//     let end = checks as i32 - 1;
//     let start: i32 = 0;
//
//     println!("score {:?}", result.score());
//     println!("query end {:?}", result.query_end());
//     println!("ref end {:?}", result.ref_end());
//     println!("query start {:?}", result.query_start());
//     println!("ref start {:?}", result.ref_start());
//     println!("cigar {:?}", result.cigar());
//     println!("cigar length {:?}", result.cigar_len());
//
//     assert_eq!(result.score(), checks);
//     assert_eq!(result.query_end(), end);
//     assert_eq!(result.ref_end(), end);
//     assert_eq!(result.query_start(), start);
//     assert_eq!(result.ref_start(), start);
//
//     Ok(())
// }
