# Contributing

Thanks for your interest in parasail-rs. There are many ways to contribute, including:

- Reporting issues
- Updating docs
- Bug fixes or feature enhancements

## Setup

The provided Nix flake can be used to set up an environment with the necessary development dependencies. To get started,
install the Nix package manager (The [Determinate Nix installer](https://determinate.systems/nix-installer/) is recommended,
but the default [Nix installer](https://nixos.org/download/#nix-install-linux) should also work fine). An environment
can be created by running `nix develop` from the project root.

Alternatively, you can set up build dependencies manually. [Rust](https://www.rust-lang.org/) (and a C compiler) are
required to build parasail-rs. CMake >= 3.5 is also required if the parasail library hasn't been compiled already.

On Ubuntu (and other Debian-based distributions), you can install a C compiler and CMake with:

```bash
sudo apt install build-essential cmake
```

## Making Changes

1. Set up a local version of parasail-rs
    - Make sure any dev dependencies are installed (if needed for your change)
    - Fork the parasail-rs repository
    - Clone the forked repository to a local directory
2. Checkout branch that should be used for the change (most likely this will be `main`)
3. Create and switch to a new Git branch
   - Ideally the branch name hints at the change you're making (e.g., `update-alignment-docs`)
4. Make and commit any changes to the local repository
   - We recommend using clippy to catch common mistakes or easy fixes
5. Create a pull request with a brief description of the proposed changes
   - Ideally using the local branch name as the base branch

## Local Testing

We encourage adding any necessary tests with your pull request.
Tests can be run locally with `cargo test [OPTIONS] [TESTNAME] ][-- [ARGS]...]` as usual.

## Documentation

To preview changes to documentation locally, use `cargo-doc`. For example,

```bash
cargo doc --open # builds and opens the docs in the default browser
```

