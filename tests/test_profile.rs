extern crate parasail_rs;

use std::ffi::CString;

use libparasail_sys::parasail_lookup_pcreator;
use parasail_rs::{Error, InstructionSet, Matrix, Profile};

const QUERY: &[u8] = b"ATGTAGTAA";

#[test]
fn default_profile_creation() -> Result<(), Error> {
    let matrix = Matrix::default();
    Profile::new(QUERY, &matrix).build()?;
    Ok(())
}

#[test]
fn profile_creation_with_stats() -> Result<(), Error> {
    let matrix = Matrix::default();
    Profile::new(QUERY, &matrix).use_stats().build()?;
    Ok(())
}

#[test]
fn profile_lookup() -> Result<(), Error> {
    // parasail_profile_create_sat()
    let fn_name = CString::new("parasail_profile_create_sat").unwrap();
    // let fn_name_ptr = fn_name.as_ptr();
    let profile_fn = unsafe { parasail_lookup_pcreator(fn_name.as_ptr()) };
    println!("{:?}", profile_fn);
    Ok(())
}

#[test]
fn profile_creation_with_solution_width() -> Result<(), Error> {
    let matrix = Matrix::default();
    let solution_widths: Vec<u8> = vec![8, 16, 32, 64];
    solution_widths.iter().for_each(|width| {
        let profile = Profile::new(QUERY, &matrix).solution_width(*width).build();
        assert!(profile.is_ok())
    });
    Ok(())
}

#[cfg(target_arch = "x86_64")]
#[test]
fn profile_x86_instructions() -> Result<(), Error> {
    let matrix = Matrix::default();
    Profile::new(QUERY, &matrix)
        .instruction_set(InstructionSet::SSE2)
        .build()?;
    Profile::new(QUERY, &matrix)
        .instruction_set(InstructionSet::SSE41)
        .build()?;
    Profile::new(QUERY, &matrix)
        .instruction_set(InstructionSet::AVX2)
        .build()?;
    Ok(())
}

#[cfg(target_arch = "aarch64")]
#[test]
fn profile_aarch64_instructions() -> Result<(), Error> {
    let matrix = Matrix::default();
    if std::arch::is_aarch64_feature_detected!("neon") {
        Profile::new(QUERY, &matrix)
            .instruction_set(InstructionSet::Neon)
            .build()?;
    }
    Ok(())
}

#[cfg(target_arch = "powerpc")]
#[test]
fn profile_powerpc_instructions() -> Result<(), Error> {
    let matrix = Matrix::default();
    if std::arch::is_powerpc64_feature_detected!(""altivec"") {
        Profile::new(QUERY, &matrix)
            .instruction_set(InstructionSet::Altivec)
            .build()?;
    }

    Ok(())
}
