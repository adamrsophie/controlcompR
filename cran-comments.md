## Submission

This is the first submission of controlcompR.

## Test environments

* GitHub Actions with `R CMD check --as-cran`: macOS (R release),
  Windows (R release), Ubuntu (R release and R devel).
* GitHub Actions, Ubuntu R devel, with the PDF and HTML manual, the CRAN
  incoming feasibility check (remote checks and aspell enabled), spelling
  and URL checks.
* GitHub Actions, Ubuntu R release, with only the hard dependencies
  installed (`_R_CHECK_FORCE_SUGGESTS_=false`).
* R-hub v2 (R devel): linux, windows, macos-arm64, valgrind, clang-asan,
  clang-ubsan and atlas.
* win-builder: R devel and R release 4.6.2 (1 note each, as below).

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.
* Possibly misspelled words in DESCRIPTION: these are author names in the
  method references (Bruinsma, Steinbuch, Stellato, Van Dooren, Van
  Overschee, De Moor, et al.) and standard control theory terms and
  acronyms (Kalman, Luenberger, Lyapunov, Riccati, Schur, Tustin,
  bumpless, controllability, observability, ARX, LQG, LQR, MPC, PID). All acronyms are
  spelled out in the Description.

## Notes for reviewers

* The package contains C code that calls LAPACK (`dgees`, `dgges`,
  `dtrsyl`) for ordered real Schur and generalized Schur decompositions,
  which base R does not expose with eigenvalue reordering.
* Under the clang sanitizer builds, function pointer type warnings are reported from the
  C sources of the 'osqp' dependency (`osqp_api.c`, `auxil.c`,
  `polish.c`), not from this package's code.
