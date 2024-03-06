# Change Log

All notable changes will be documented here in reverse chronological order the headers \<VERSION\> - <YY.MM.DD>.

## 0.4.0 - 2024.03.xx

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

--

TODO:
Implement AlignResult methods:
    - [ ] get_trace_ins_table
    - [ ] get_trace_del_table
- [ ] convert panics to better errors
    - [ ] use results for handling errors
- [ ] convert unwraps for better error handling
- [x] using new builder pattern for Aligner
- [x] global, local, and semi_global alignment stats functions
- [x] implement different semi_global variants
- [x] local alignment implementation
- [x] semi-global alignment implementation
- [x] tests for scoring matrices and alignments
- [x] query profiles and corresponding functions
- [x] Work on adding safety to bindings

