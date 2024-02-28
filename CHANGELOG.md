# Change Log

## 0.2.1 - 2024.02.28

- Use Arc for profile and matrix instead of Rc for thread safety

## 0.2.0 - 2024.02.25

- Use new builder pattern for aligner.
- Better Matrix and Profile handling with Rc.
- Methods for AlignResults struct.

## 0.1.0 - 2024.02.21

Initial release.

- Implements global, semi-global, and local alignment.
- Supports query profiles or one-off alignments.

--

TODO:
- [x] using new builder pattern for Aligner
- [ ] global, local, and semi_global alignment stats functions
- [x] implement different semi_global variants
- [x] local alignment implementation
- [x] semi-global alignment implementation
- [ ] tests for scoring matrices and alignments
- [x] query profiles and corresponding functions
- [x] Work on adding safety to bindings

