# Contributing to controlcompR

Bug reports and pull requests are welcome at
<https://github.com/adamrsophie/controlcompR>.

- Report numerical problems with a minimal reproducible example, the
  expected result and, when possible, the result of another tool
  (python-control, MATLAB, ControlSystems.jl).
- New algorithms need tests that check mathematical properties
  (residuals, invariances, agreement with an independent method), not
  only stored numbers.
- Functions must refuse problems they cannot solve reliably, with an
  informative error.
- Use American English and roxygen2 documentation for all exported
  functions.
