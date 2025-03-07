extern crate parasail_rs;

use parasail_rs::{Aligner, Matrix, Profile};

#[test]
fn build_aligner() {
    assert!(Aligner::new().build().is_ok()); // default is global aligner
    assert!(Aligner::new().local().build().is_ok());
    assert!(Aligner::new().semi_global().build().is_ok());
}

#[test]
fn default_alignment() -> Result<(), Box<dyn std::error::Error>> {
    let query = b"ACTTAGCTT";
    let reference = b"ACTTAGCTT";
    let global_aligner = Aligner::new().build()?;
    let local_aligner = Aligner::new().local().build()?;
    let semi_global_aligner = Aligner::new().semi_global().build()?;

    let global_result = global_aligner.align(query, reference)?;
    let local_result = local_aligner.align(query, reference)?;
    let semi_global_result = semi_global_aligner.align(query, reference)?;

    assert_eq!(local_result.score()?, query.len() as i32);
    assert_eq!(global_result.score()?, query.len() as i32);
    assert_eq!(semi_global_result.score()?, query.len() as i32);

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

    assert_eq!(global_result.score()?, query.len() as i32);
    assert_eq!(semi_global_result.score()?, query.len() as i32);
    assert_eq!(local_result.score()?, query.len() as i32);

    Ok(())
}
