# testthat is a suggested package: run the tests only when it is available
if (requireNamespace("testthat", quietly = TRUE)) {
  library(testthat)
  library(controlcompR)
  test_check("controlcompR")
}
