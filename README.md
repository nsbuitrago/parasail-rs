# parasail-rs

[![GitHub Actions Workflow
Status](https://img.shields.io/github/actions/workflow/status/nsbuitrago/parasail-rs/test.yml)](https://github.com/nsbuitrago/parasail-rs/actions)
[![docs.rs](https://img.shields.io/docsrs/parasail-rs)](https://docs.rs/parasail-rs/latest/parasail-rs)
[![Crates.io
Version](https://img.shields.io/crates/v/parasail-rs)](https://crates.io/crates/parasail-rs)

This crate provides Rust wrapper for
[parasail](https://github.com/jeffdaily/parasail/tree/master), a SIMD pairwise
sequence alignment C library. Note that this crate is still under development
and is unstable.

## Usage

### Installation

Run the following Cargo command in your project directory:

```bash
cargo add parasail_rs
```

Note that parasail-rs depends on libparasail-sys which will either use an
already installed system parasail library or build from source. For more
information, please see
[libparasail-sys](https://github.com/nsbuitrago/libparasail-sys).

### Examples

#### Basic usage

For one-off alignments:

```rust
use parasail_rs::Aligner;

let query = b"ACGT";
let reference = b"TAGACGTTTA";
let aligner = Aligner::new().build();

aligner.align(Some(query), reference)?;
```

#### Creating scoring matrices

Using query profile:

```rust
use parasail_rs::{Matrix, Aligner, Profile};

// ...

let query = b"ACGT";
let ref_1 = b"ACGTAACGTACA";
let ref_2 = b"TGGCAAGGTAGA";

let use_stats = true;
let query_profile = Profile::new(query, use_stats, &Matrix::default())?;
let aligner = Aligner::new()
    .profile(query_profile)
    .build();

let result_1 = aligner.align(None, ref_1)?;
let result_2 = aligner.align(None, ref_2)?;

println!("Score 1: {}", result_1.get_score());
println!("Score 2: {}", result_2.get_score());
```

## Contributing

Contributions are more than welcome. Please open an issue for any feedback or questions.

## Citations

If needed, please cite the following paper:

[Daily, Jeff. (2016). Parasail: SIMD C library for global, semi-global, and local pairwise sequence alignments. BMC Bioinformatics, 17(1), 1-11. doi:10.1186/s12859-016-0930-z](https://doi.org/10.1186/s12859-016-0930-z)

## License

parasail-rs and libparasail-sys are licensed under the BSD-3-Clause license. The
original parasail C library is licensed under a similar Battelle style BSD
license.
