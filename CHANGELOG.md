# Change Log

All notable changes will be documented here in reverse chronological order the headers \<VERSION\> - <YY.MM.DD>.

## 0.9.0 - 2026.01.04

### Added

- New `Table` struct providing access to alignment tables (scores, similar, length).
  - `.get(row, col)` for bounds-checked access to individual cells.
  - `.rows()` and `.cols()` for table dimensions.
  - `.last()` for the last cell (typically the final alignment score).
  - `.as_slice()` for raw data access.
  - `Display` implementation for pretty-printing tables.
- `TracebackTable` struct providing access to traceback table.
  - `.get(row, col)` for bounds-checked access to individual cells returning simple flags
    including 'DIAG', 'INS', and 'DEL'.
  - `.get_detailed(row, col)` for bounds-checked access to individual cells returning all
    parasail flags (such as 'INS_E`, 'DEL_F', etc.)
  - `.rows()` and `.cols()` for table dimensions.
  - `.as_slice()` for raw data access.
  - `fmt::Display` implementation for pretty-printing tables with simple flags, and
    `fmt::Debug` implementation for printing detailed parasail flags. 
- `bitflags` dependency for ergonmic TraceFlags generation.

### Breaking Changes

- Fixed table getter methods to return full tables instead of raw pointer:
  - `get_score_table()` now returns `Result<Table>` instead of `Result<i32>`
  - `get_matches_table()` now returns `Result<Table>` instead of `Result<i32>`
  - `get_similar_table()` now returns `Result<Table>` instead of `Result<i32>`
  - `get_length_table()` now returns `Result<Table>` instead of `Result<i32>`
  - `get_trace_table()` now returns `Result<TracebackTable>` instead of `Result<i32>`

### Fix

- Properly expose `AlignerBuilder`, `Alignment`, and `Traceback` structs.

### Update

- bump deps:
  - thiserror 2.0.12 -> 2.0.17
  - derive_more 2.0.1 -> 2.1.1
  - libc 0.2.174 -> 0.2.179
  - log 0.4.27 -> 0.4.29
  - libparasail-sys 0.1.12 -> 0.2.0
- parasail-rs error enums moved to use derive_more.
- Remove thiserror dep.

## 0.8.1 - 2025.08.01

### Fix

- ARM Linux architecture support

## 0.8.0 - 2025.06.24

### Changes

- Add central Error enum and Result type
- Add SolutionWidth and InstructionSet enums for building specific profiles
- Add ProfileBuilder struct for more control over profile configuration
- Rename AlignResult -> Alignment
  - Add query and ref length to the Alignment struct

### Fix

- output row/col slice for row/col alignment methods (thanks to @ctsa)

### Update

- bump deps:
  - libparasail-sys 0.1.10 -> 0.1.11
  - libc 0.2.172 -> 0.2.174
  - thiserror 1.0.69 -> 2.0.12
- multithreaded alignment test

## 0.7.8 - 2025.06.03

### Update

- bump deps. libparasail-sys v0.1.10 removes support for cmake < 3.5
    - libparasail-sys 0.1.9 -> 0.1.10
    - libc 0.2.169 -> 0.2.172
    - log 0.4.25 -> 0.4.27

## 0.7.7 - 2025.01.19

### Update

- bump dependencies. libparasail-sys v0.1.9 now skips building tests and the
aligner application by default.

## 0.7.6 - 2024.11.18

### Fix

- bindgen bump issue with unstable library feature. (see
[#2](https://github.com/nsbuitrago/libparasail-sys/pull/2)).

## 0.7.5 - 2024.11.1

### Update

- bump deps

## 0.7.4 - 2024.08.24

### Bug Fixes

- Add `solution_width` method to `AlignerBuilder` (#8)
- Add copy trait to aligner (#7)

## 0.7.3 - 2024.04.22

### Features

- Add SSW library emulation. The `ssw` method on `Aligner` performs
striped Smith-Waterman local alignment using SSE2 instructions. Support for
profiles will be supported in a later release.
- A new `SSWResult` has been added and is returned from successful `ssw`
alignment. With SSWResult, you can get the primary alignment score, along with
start and end locations of the alignment on the query or reference.
- Add banded global (Needleman-Wuncsh) alignment. The `banded_ssw` method
on `Aligner` performs a banded global alignment.
- Add `bandwith` method on `AlignerBuilder` to set the band with for banded nw
alignment.

### Improvements

- Updated doc strings in lib.rs

## 0.7.2 - 2024.04.17

### Bug Fixes

- AlignerFn enum uses *const c_char in expected fn signature instead of *const i8 as before
- With fix above, we can use as_ptr() as in the query/ref and parasail fn name as
before without casting to *const i8!

## 0.7.1 - 2024.04.17

- Upgrade dependencies. Notably libparasail-sys -> 0.1.6

### Bug Fixes

- Use *const i8 casting for passing fn_name into `parasail_lookup_function`
and query/reference to alignment fn. This attempts to solve a compilation
issue on aarch64 linux with gcc-multilib available.

## 0.7.0 - 2024.03.20

### Features

- Add Display trait on Matrix for easy printing to stdout.

### Breaking Changes

- Matrix methods returning results:
    - `create`, `from`, and `create_pssm`, `from_file` now return `Result<Self, MatrixError>`.
- Profile `create` method now returns `Result<Self, ProfileError>`.
- `align` method returns `Result<AlignResult, AlignError>`.
- `convert_square_to_pssm` renamed to `to_pssm` for succinctness.

### Bug Fixes

- `set_value` method can only modify valid matrix indices.

## 0.6.0 - 2024.03.18

### Features

- Add `set_value` method to change values in substitution matrices.
- Add methods on AlignResult:
    - `print_traceback` (print alignment strings and statistics).
    - `get_traceback_strings` (return a Traceback struct that includes query, ref, and comparison strings).
    - `get_cigar` (get CIGAR strings).

## Breaking Changes

- Removed `mode` method from AlignerBuilder. Use `global`, `semi_global`, or `local` methods instead.
- Removed `vec_strategy` method from AlignerBuilder. Use `striped`, `scan`, and `diag` methods instead.

## 0.5.1 - 2024.03.13

### Features

- Add `global`, `semi_global`, and `local` methods to `AlignerBuilder` to set alignment algorithm.
- Add `striped`, `scan`, and `diag` methods to `AlignerBuilder` to set vectorization method.

## 0.5.0 - 2024.03.11

### Breaking Changes

- Aligner builder takes now has a mode method to set alignment algorithm (nw, sg, sw).
- Many wrapper alignment methods removed in favor of single align method.

### Bug Fixes

- constructing parasail fn names from Aligner builder to avoid building on every alignment.

## 0.4.0 - 2024.03.06

### Features

- Better support for allowing various gap configurations in semi-global alignments [\#1].

### Breaking Changes

- removed allow_gaps
- use allow_query_gaps and allow_ref_gaps in AlignerBuilder / Aligner for semi-global alignment.

### Bug Fixes

- use_table and use_trace disable each other if called to be mutually exclusive [\#2].
- use_table and use_stats disable each other if called to be mutually exclusive [\#2].

## 0.3.0 - 2024.03.03

- Note that currenlty use_trace and use_stats can both be set to true. However, traceback is only compatible without stats.

### Breaking Changes

- use_stats method changed to is_stats on AlignResult.
- use_stats and use_table methods on AlignerBuilder do not take a bool. Simply
call the methods to set them to true. By default, these extra params for alignment
are not used.

### Features

New methods for AlignResults:

- [x] get_score_table
- [x] get_matches_table
- [x] get_similar_table
- [x] get_length_table
- [x] get_score_row
- [x] get_matches_row
- [x] get_similar_row
- [x] get_length_row
- [x] get_score_col
- [x] get_matches_col
- [x] get_similar_col
- [x] get_length_col
- [x] get_trace_table
- [x] is_blocked
- [x] is_stats_table
- [x] is_stats_rowcol
- [x] is_rowcol
- [x] is_stats_trace

### Bug Fixes

- Check that aligner was initialized with use_stats = true before calling methods that requires additional stats:
    - get_matches
    - get_length
- Check that use_stats for profile and aligner match before aligning
- Avoid calling `parasail_matrix_free` on drop for builtin matrices

## 0.2.3 - 2024.02.28

### Changes

- Pass use_stats to aligner builder to use alignment functions that return additional statistics

### Bug Fixes

- Some AlignResult getter functions returning wrong boolean

## 0.2.2 - 2024.02.28

### Bug Fixes

- Add Send and Sync traits to profile and matrix structs

## 0.2.1 - 2024.02.28

### Changes

- Use Arc for profile and matrix instead of Rc for thread safety

## 0.2.0 - 2024.02.25

### Changes

- Use new builder pattern for aligner.
- Better Matrix and Profile handling with Rc.
- Methods for AlignResults struct.

## 0.1.0 - 2024.02.21

### Initial release.

- Implements global, semi-global, and local alignment.
- Supports query profiles or one-off alignments.
