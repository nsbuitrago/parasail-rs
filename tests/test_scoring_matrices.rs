use parasail_rs::{ScoringMatrix, Aligner};
use std::io;

#[test]
pub fn test_scoring_matrix() -> Result<(), io::Error> {

    let sm = ScoringMatrix::from("blosum62");
    unsafe {
        println!("{:?}", (*sm.scores));
    }

    let query = b"ACGT";
    let target = b"ACGT";
    let aligner = Aligner::new(sm, 5, 2);

    let result = aligner.global(query, target, "striped");
    println!("{:?}", result);

    let score = result.get_score();
    println!("Score: {}", score);

    Ok(())
}
