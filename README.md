# parasail-rs

![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/nsbuitrago/parasail-rs/test.yml) [![docs.rs](https://img.shields.io/docsrs/parasail-rs)](https://docs.rs/parasail-rs/latest/parasail_rs/index.html) [![Crates.io Version](https://img.shields.io/crates/v/parasail-rs)](https://crates.io/crates/parasail-rs) [![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/nsbuitrago/parasail-rs)

This crate provides safe Rust bindings and a wrapper to [parasail](https://github.com/jeffdaily/parasail/tree/master), a SIMD pairwise sequence alignment C library. Note that this crate is still under development and is unstable.

## Usage

### Installation

Run the following Cargo command in your project directory:

```bash
cargo add parasail-rs
```

Note that parasail-rs depends on libparasail-sys which will either use an already installed system parasail library or build from source. For more information, please see [libparasail-sys](https://github.com/nsbuitrago/libparasail-sys).

### Examples

#### Basic usage:

For one-off alignments:

```rust
use parasail_rs::prelude::Aligner;

// ...

let query = b"ACGT";
let reference = b"ACGT";
let aligner = Aligner::new().build();

aligner.align(Some(query), reference)?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

When using striped or scan vectorization strategies, some performance may
be gained by reusing the query sequence. This can be done by creating a
query profile and reusing it for multiple alignments.

```rust
use parasail_rs::prelude::{Matrix, Aligner, Profile};

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
# Ok::<(), Box<dyn std::error::Error>>(())
```

## Contributing

Contributions are more than welcome. Please open an issue for any feedback or questions. See the [contributing guide](CONTRIBUTING.md)
for more details.

## Citations

If needed, please cite the following paper:

[Daily, Jeff. (2016). Parasail: SIMD C library for global, semi-global, and local pairwise sequence alignments. BMC Bioinformatics, 17(1), 1-11. doi:10.1186/s12859-016-0930-z](https://doi.org/10.1186/s12859-016-0930-z)

## License

parasail-rs and libparasail-sys are licensed under the BSD-3-Clause license. The original parasail C library is licensed under a similar Battelle style BSD license.
