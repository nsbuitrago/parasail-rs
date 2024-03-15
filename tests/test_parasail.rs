use parasail_rs::{Aligner, Matrix, Profile};
use std::thread;

#[test]
pub fn matrix_construction() {
    // default matrix
    Matrix::default();

    // custom matrix
    Matrix::create(b"ACGT", 3, -2);

    // from name
    Matrix::from("blosum62");

    // square matrix from file
    Matrix::from_file("./tests/square.txt");

    // pssm from file
    Matrix::from_file("./tests/pssm.txt");

    // PSSM
    // let values = vec![1, 2, 3, 4];
    // let rows = vec![b'A', b'C'];
    // Matrix::create_pssm("AC", values, rows);

    // convert square to pssm
}

#[test]
pub fn profile_construction() {
    let query = b"ATGGCACTATAA";
    Profile::new(query, false, &Matrix::default());

    // with stats
    Profile::new(query, true, &Matrix::default());
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
        .vec_strategy("striped")
        .use_stats()
        .build();
}

#[test]
pub fn global_alignment() {
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().build();
    let result = aligner.align(Some(query), reference);

    let checks = query.len() as i32;

    assert_eq!(result.get_score(), checks);
    assert_eq!(result.get_end_query(), checks - 1);
    assert_eq!(result.get_end_ref(), checks - 1);
    assert!(result.is_global());
    assert!(!result.is_local());
    assert!(!result.is_semi_global());
    assert!(result.is_striped());
}

#[test]
pub fn semi_global_alignment() {
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().semi_global().build();
    let result = aligner.align(Some(query), reference);

    let checks = query.len() as i32;

    assert_eq!(result.get_score(), checks);
    assert_eq!(result.get_end_query(), checks - 1);
    assert_eq!(result.get_end_ref(), checks - 1);
    assert!(!result.is_global());
    assert!(!result.is_local());
    assert!(result.is_semi_global());
    assert!(result.is_striped());
}

#[test]
pub fn local_alignment() {
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().local().build();
    let result = aligner.align(Some(query), reference);

    let checks = query.len() as i32;

    assert_eq!(result.get_score(), checks);
    assert_eq!(result.get_end_query(), checks - 1);
    assert_eq!(result.get_end_ref(), checks - 1);
    assert!(!result.is_global());
    assert!(result.is_local());
    assert!(!result.is_semi_global());
    assert!(result.is_striped());
}

#[test]
pub fn global_with_stats() -> Result<(), Box<dyn std::error::Error>> {
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().use_stats().build();
    let result = aligner.align(Some(query), reference);

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
    let aligner = Aligner::new().semi_global().use_stats().build();
    let result = aligner.align(Some(query), reference);

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
    let aligner = Aligner::new().local().use_stats().build();
    let result = aligner.align(Some(query), reference);

    let checks = query.len() as i32;
    let n_matches = result.get_matches()?;
    let align_len = result.get_length()?;

    assert_eq!(n_matches, checks);
    assert_eq!(align_len, checks);

    Ok(())
}

#[test]
pub fn score_table() -> Result<(), Box<dyn std::error::Error>> {
    // one-off alignment wihthout stats
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().use_table().build();
    let default_score = 1;

    let result = aligner.align(Some(query), reference);
    assert!(result.is_table());
    assert!(!result.is_stats());
    assert!(!result.is_stats_table());
    assert_eq!(result.get_score_table()?, default_score);

    // one-off alignment with stats
    let aligner = Aligner::new().use_stats().use_table().build();

    let result = aligner.align(Some(query), reference);
    assert!(result.is_stats());
    assert!(result.is_stats_table());
    assert!(result.is_table());
    assert_eq!(result.get_score_table()?, default_score);

    // alignment with profile, without stats
    let custom_score = 3;
    let matrix = Matrix::create(b"ACGT", custom_score, -2);
    let profile = Profile::new(query, false, &matrix);

    let aligner_w_profile = Aligner::new().profile(profile).use_table().build();

    let result_w_profile = aligner_w_profile.align(None, reference);

    assert!(result_w_profile.is_table());
    assert!(!result_w_profile.is_stats());
    assert!(!result_w_profile.is_stats_table());
    assert_eq!(result_w_profile.get_score_table()?, custom_score);

    // // alignment with profile, with stats
    let profile = Profile::new(query, true, &matrix);
    let aligner_w_profile = Aligner::new()
        .profile(profile)
        .use_stats()
        .use_table()
        .build();

    let result_w_profile = aligner_w_profile.align(None, reference);

    assert!(result_w_profile.is_stats());
    assert!(result_w_profile.is_stats_table());
    assert!(result_w_profile.is_table());
    assert_eq!(result_w_profile.get_score_table()?, custom_score);

    Ok(())
}

#[test]
pub fn matches_table() -> Result<(), Box<dyn std::error::Error>> {
    // one-off alignment wihthout stats
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().use_table().use_stats().build();
    let default_score = 1;

    let result = aligner.align(Some(query), reference);
    assert!(result.is_table());
    assert!(result.is_stats());
    assert!(result.is_stats_table());
    assert_eq!(result.get_matches_table()?, default_score);

    Ok(())
}

#[test]
pub fn similar_table() -> Result<(), Box<dyn std::error::Error>> {
    // one-off alignment wihthout stats
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().use_table().use_stats().build();

    let result = aligner.align(Some(query), reference);
    assert!(result.is_table());
    assert!(result.is_stats());
    assert!(result.is_stats_table());
    println!("similar table: {:?}", result.get_similar_table()?);

    Ok(())
}

#[test]
pub fn length_table() -> Result<(), Box<dyn std::error::Error>> {
    // one-off alignment wihthout stats
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().use_table().use_stats().build();

    let result = aligner.align(Some(query), reference);
    assert!(result.is_table());
    assert!(result.is_stats());
    assert!(result.is_stats_table());
    println!("Length table: {:?}", result.get_length_table()?);

    Ok(())
}

#[test]
pub fn score_row() -> Result<(), Box<dyn std::error::Error>> {
    // one-off alignment wihthout stats
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().use_last_rowcol().use_stats().build();

    let result = aligner.align(Some(query), reference);
    assert!(result.is_stats_rowcol());
    assert!(result.is_stats());
    assert!(!result.is_stats_table());
    println!("Score row: {:?}", result.get_score_row()?);

    Ok(())
}

#[test]
pub fn matches_row() -> Result<(), Box<dyn std::error::Error>> {
    // one-off alignment wihthout stats
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().use_last_rowcol().use_stats().build();

    let result = aligner.align(Some(query), reference);
    assert!(result.is_stats_rowcol());
    assert!(result.is_stats());
    assert!(!result.is_stats_table());
    println!("Matches row: {:?}", result.get_matches_row()?);

    Ok(())
}

#[test]
pub fn similar_row() -> Result<(), Box<dyn std::error::Error>> {
    // one-off alignment wihthout stats
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().use_last_rowcol().use_stats().build();

    let result = aligner.align(Some(query), reference);
    assert!(result.is_stats_rowcol());
    assert!(result.is_stats());
    assert!(!result.is_stats_table());
    println!("Similar row: {:?}", result.get_similar_row()?);

    Ok(())
}

#[test]
pub fn length_row() -> Result<(), Box<dyn std::error::Error>> {
    // one-off alignment wihthout stats
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().use_last_rowcol().use_stats().build();

    let result = aligner.align(Some(query), reference);
    assert!(result.is_stats_rowcol());
    assert!(result.is_stats());
    assert!(!result.is_stats_table());
    println!("Length row: {:?}", result.get_length_row());

    Ok(())
}

#[test]
pub fn score_col() -> Result<(), Box<dyn std::error::Error>> {
    // one-off alignment wihthout stats
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().use_last_rowcol().use_stats().build();

    let result = aligner.align(Some(query), reference);
    assert!(result.is_stats_rowcol());
    assert!(result.is_stats());
    assert!(!result.is_stats_table());
    println!("Score col: {:?}", result.get_score_col()?);

    Ok(())
}

#[test]
pub fn match_col() -> Result<(), Box<dyn std::error::Error>> {
    // one-off alignment wihthout stats
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().use_last_rowcol().use_stats().build();

    let result = aligner.align(Some(query), reference);
    assert!(result.is_stats_rowcol());
    assert!(result.is_stats());
    assert!(!result.is_stats_table());
    println!("Matches col: {:?}", result.get_matches_col()?);

    Ok(())
}

#[test]
pub fn similar_col() -> Result<(), Box<dyn std::error::Error>> {
    // one-off alignment wihthout stats
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().use_last_rowcol().use_stats().build();

    let result = aligner.align(Some(query), reference);
    assert!(result.is_stats_rowcol());
    assert!(result.is_stats());
    assert!(!result.is_stats_table());
    println!("Similar col: {:?}", result.get_similar_col()?);

    Ok(())
}

#[test]
pub fn length_col() -> Result<(), Box<dyn std::error::Error>> {
    // one-off alignment wihthout stats
    let query = b"ACGT";
    let reference = b"ACGT";
    let aligner = Aligner::new().use_last_rowcol().use_stats().build();

    let result = aligner.align(Some(query), reference);
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
    let aligner = Aligner::new().use_trace().build();
    let result = aligner.align(Some(query), reference);
    assert!(result.is_trace());
    println!("Trace table: {:?}", result.get_trace_table()?);

    Ok(())
}

// #[test]
// pub fn trace_ins_table() -> Result<(), Box<dyn std::error::Error>> {
//     let query = b"ACGTA";
//     let reference = b"ACGTTA";
//     let aligner = Aligner::new().use_trace().build();
//     let result = aligner.global(query, reference);
//     assert!(result.is_trace());
//     // println!("Trace ins table: {:?}", result.get_trace_ins_table()?);
//
//     Ok(())
// }

// #[test]
// pub fn trace_del_table() -> Result<(), Box<dyn std::error::Error>> {
//     let query = b"ACGT";
//     let reference = b"ACGT";
//     let aligner = Aligner::new().use_trace().build();
//     let result = aligner.global(query, reference);
//     assert!(result.is_trace());
//     println!("Trace del table: {:?}", result.get_trace_del_table()?);
//
//     Ok(())
// }

#[test]
pub fn global_with_profile() {
    let query = b"ACGT";
    let reference = b"ACGT";
    let matrix = Matrix::default();
    let profile = Profile::new(query, true, &matrix);

    let aligner = Aligner::new()
        .profile(profile)
        .use_stats()
        .vec_strategy("striped")
        .build();

    let result = aligner.align(None, reference);
    assert!(result.is_global());
    assert!(result.is_striped());
    assert!(result.is_stats());
    assert!(!result.is_local());
    assert!(!result.is_semi_global());
}

#[test]
pub fn semi_global_with_profile() {
    let query = b"ACGT";
    let reference = b"ACGT";
    let matrix = Matrix::default();
    let profile = Profile::new(query, true, &matrix);

    let aligner = Aligner::new()
        .profile(profile)
        .use_stats()
        .vec_strategy("striped")
        .semi_global()
        .build();

    let result = aligner.align(None, reference);
    assert!(result.is_semi_global());
    assert!(result.is_striped());
    assert!(result.is_stats());
    assert!(!result.is_local());
    assert!(!result.is_global());
}

#[test]
pub fn local_with_profile() {
    let query = b"ACGT";
    let reference = b"ACGT";
    let matrix = Matrix::default();
    let profile = Profile::new(query, true, &matrix);

    let aligner = Aligner::new()
        .profile(profile)
        .use_stats()
        .vec_strategy("striped")
        .local()
        .build();

    let result = aligner.align(None, reference);
    assert!(result.is_local());
    assert!(result.is_striped());
    assert!(result.is_stats());
    assert!(!result.is_global());
    assert!(!result.is_semi_global());
}

#[test]
pub fn multithread_global_alignment() {
    let query = b"ACGT";
    let refs = vec![b"ACGT", b"ACGT"];
    let matrix = Matrix::default();
    let profile = Profile::new(query, true, &matrix);

    let aligner = Aligner::new().profile(profile).use_stats().build();

    thread::spawn(move || {
        for reference in refs {
            let result = &aligner.align(None, reference);
            let score = result.get_score();
            assert_eq!(score, query.len() as i32);
        }
    })
    .join()
    .unwrap();
}
