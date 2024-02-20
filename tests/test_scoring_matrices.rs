use parasail_rs::{ScoringMatrix, Aligner, Profile};
use std::io;

#[test]
pub fn test_scoring_matrix() -> Result<(), io::Error> {

    let sm = ScoringMatrix::from("blosum62");
    unsafe {
        println!("{:?}", (*sm.scores));
    }

    let query = b"ACGT";
    let target = b"ACGT";
    let vector_strategy = String::from("striped");
    // let aligner = Aligner::new(sm, 5, 2, vector_strategy);
    let profile = Profile::new(query, false, sm);
    let aligner = Aligner::with_profile(profile, 5, 2, vector_strategy);

    let result = aligner.global_with_profile(target);
    println!("{:?}", result);

    let score = result.get_score();
    println!("Score: {}", score);

    Ok(())
}
