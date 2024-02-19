# parasail-rs

This crate provides safe Rust bindings to [parasail](), a vectorized pairwise sequence alignment C library.

## Usage

### Installation

Add the latest version of `parasail-rs` to your `Cargo.toml` by running `cargo add parasail-rs`. Or, add directory to `Cargo.toml`:

```toml

[dependencies]
parasail-rs = "0.1.0"

```

Note that parasail-rs depends on libparasail-sys which will either use an already installed system parasail library or build from source. For more information, please see [libparasail-sys](https://gitlab.com/nsbuitrago/libparasail-sys).

### Examples

For the moment, there are some examples in the doc string of lib.rs. These will be documented further when the crate is published on crates.io

## Contributing

Contributions are more than welcome. For other feedback or suggestions, please open an issue or send an email to nsb5 [at] rice.edu.

## License

parasail-rs and libparasail-sys are licensed under the BSD-3-Clause license. The original parasail C library is licensed under a similar Battelle style BSD license.

Nicolas S. Buitrago \<nsb5 [at] rice.edu\>

