extern crate parasail_rs;
use parasail_rs::{Aligner, Matrix};

#[test]
fn global_alignment() -> Result<(), Box<dyn std::error::Error>> {
    let query = b"ACTTAGCTT";
    let reference = b"ACTTAGCTT";
    let matrix = Matrix::default();
    let aligner = Aligner::new().build()?;

    let result = aligner.align(query, reference)?;
    //assert_eq!(result.score(), query.len());
    Ok(())
}

// fn local_alignment() -> Result<(), Box<dyn std::error::Error>> {
//     let query = b"ACTTAGCTT";
//     let reference = b"ACTTAGCTT";
//     let matrix = Matrix::default();
//     let aligner = Aligner::new().local().build()?;

//     let result = aligner.align(query, reference)?;
//     // assert_eq!(result.score(), query.len());
//     //
//     Ok(())
// }
