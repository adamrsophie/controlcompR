## Submission

This is the first submission of controlcompR.

## Test environments

* GitHub Actions: macOS (R release), Windows (R release), Ubuntu (R release
  and R devel), all with `R CMD check --as-cran`.
* GitHub Actions, Ubuntu R devel with the PDF manual and the CRAN incoming
  feasibility checks enabled.
* R-hub v2: linux, windows, macos-arm64, clang-asan, valgrind and nosuggests.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Notes for reviewers

* The package contains C code that calls LAPACK (`dgees`, `dgges`,
  `dtrsyl`) for ordered real Schur and generalized Schur decompositions,
  which base R does not expose with eigenvalue reordering.
* Words flagged as possibly misspelled in the DESCRIPTION (for example
  Riccati, Lyapunov, Tustin, Luenberger) are standard control theory terms
  and author names.
