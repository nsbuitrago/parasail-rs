# parasail-rs

This crate provides safe Rust bindings and a wrapper to [parasail](), a SIMD pairwise sequence alignment C library.

## Usage

### Installation

Add the latest version of `parasail-rs` to your `Cargo.toml` by running `cargo add parasail-rs`. Or, add directory to `Cargo.toml`:

```toml

[dependencies]
parasail-rs = "0.1.0"

```

Note that parasail-rs depends on libparasail-sys which will either use an already installed system parasail library or build from source. For more information, please see [libparasail-sys](https://gitlab.com/nsbuitrago/libparasail-sys).

### Examples

Basic usage:

 ```rust
use parasail_rs::{ScoringMatrix, Aligner};
let matrix = ScoringMatrix::create("ACGT", 1, -1);
let vector_strategy = "striped".to_string();
let query = b"ACGT";
let target = b"ACGTAACGTACA";

let aligner = Aligner::new(matrix, 5, 2, vector_strategy);
let result = aligner.global(query, target); 
 ```

Using query profiles:

```rust
use parasail_rs::{ScoringMatrix, Aligner, Profile};
let matrix = ScoringMatrix::create("ACGT", 1, -1);
let vector_strategy = "striped".to_string();
let query = b"ACGT";
let ref_1 = b"ACGTAACGTACA";
let ref_2 = b"TGGCAAGGTAGA";

let use_stats = true;
let query_profile = Profile::new(query, use_stats, matrix);
let aligner = Aligner::with_profile(query_profile, 5, 2, vector_strategy);

let result_1 = aligner.global_with_profile(ref_1);
let result_2 = aligner.global_with_profile(ref_2);
```

## Contributing

Contributions are more than welcome. Please open an issue or send an email to nsb5 [at] rice.edu for any feedback or questions.

## Citations

If needed, please cite the following paper:

[Daily, Jeff. (2016). Parasail: SIMD C library for global, semi-global, and local pairwise sequence alignments. BMC Bioinformatics, 17(1), 1-11. doi:10.1186/s12859-016-0930-z](https://doi.org/10.1186/s12859-016-0930-z)

The Parasail C library was developed by Jeff Daily.

## License

parasail-rs and libparasail-sys are licensed under the BSD-3-Clause license. The original parasail C library is licensed under a similar Battelle style BSD license.

Nicolas S. Buitrago \<nsb5 [at] rice.edu\>

