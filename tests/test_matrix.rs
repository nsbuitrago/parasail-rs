extern crate parasail_rs;
use parasail_rs::Matrix;

#[test]
fn create_builtin_matrix() {
    BUILTIN_MATRICES.iter().for_each(|builtin_name| {
        assert!(Matrix::from(builtin_name).is_ok());
    });
}

#[test]
fn create_custom_matrix() {
    assert!(Matrix::create(b"ACGT", 3, -2).is_ok())
}

#[test]
fn create_matrix_from_file() {
    assert!(Matrix::from_file("tests/data/square.txt").is_ok());
    assert!(Matrix::from_file("tests/data/pssm.txt").is_ok());
}

#[test]
fn edit_matrix() -> Result<(), Box<dyn std::error::Error>> {
    // edit a custom matrix, easy
    let mut custom_dna = Matrix::default();
    custom_dna.set_value(0, 0, 3)?;

    // edit a builtin matrix. will copy the matrix first and create a custom variant
    let mut custom_blosum = Matrix::from("blosum62")?;
    custom_blosum.set_value(0, 0, 100)?;

    assert_eq!(custom_dna.get_value(0, 0)?, 3);
    assert_eq!(custom_blosum.get_value(0, 0)?, 100);

    Ok(())
}

const BUILTIN_MATRICES: &[&str; 67] = &[
    "blosum30",
    "blosum35",
    "blosum40",
    "blosum45",
    "blosum50",
    "blosum55",
    "blosum60",
    "blosum62",
    "blosum65",
    "blosum70",
    "blosum75",
    "blosum80",
    "blosum85",
    "blosum90",
    "blosum100",
    "pam10",
    "pam20",
    "pam30",
    "pam40",
    "pam50",
    "pam60",
    "pam70",
    "pam80",
    "pam90",
    "pam100",
    "pam110",
    "pam120",
    "pam130",
    "pam140",
    "pam150",
    "pam160",
    "pam170",
    "pam180",
    "pam190",
    "pam200",
    "pam210",
    "pam220",
    "pam230",
    "pam240",
    "pam250",
    "pam260",
    "pam270",
    "pam280",
    "pam290",
    "pam300",
    "pam310",
    "pam320",
    "pam330",
    "pam340",
    "pam350",
    "pam360",
    "pam370",
    "pam380",
    "pam390",
    "pam400",
    "pam410",
    "pam420",
    "pam430",
    "pam440",
    "pam450",
    "pam460",
    "pam470",
    "pam480",
    "pam490",
    "pam500",
    "dnafull",
    "nuc44",
];
