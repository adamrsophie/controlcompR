# Compare controlcompR with reference values from python-control/slycot.
# Run after reference.py from the repository root:
#   Rscript benchmarks/reference/compare.R
library(controlcompR)
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required")
prob <- jsonlite::read_json("benchmarks/reference/problems.json", simplifyVector = TRUE)
ref <- jsonlite::read_json("benchmarks/reference/reference-values.json",
                           simplifyVector = FALSE)
# JSON arrays: a list of rows becomes a row-major matrix, a list of scalars a vector
as_num <- function(x) {
  if (!is.list(x)) return(as.numeric(x))
  # complex arrays are stored as {"re": ..., "im": ...}; values compared here are real
  if (identical(sort(names(x)), c("im", "re"))) {
    im <- as_num(x$im)
    if (any(abs(im) > 1e-12)) stop("Unexpected complex reference value.")
    return(as_num(x$re))
  }
  if (all(vapply(x, is.list, TRUE))) return(do.call(rbind, lapply(x, as_num)))
  as.numeric(unlist(x))
}
m <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  x
}
rows <- list()
add <- function(name, ours, theirs) {
  theirs <- as_num(theirs)
  # vectors may arrive as 1 x n or n x 1 matrices; compare those as vectors
  if (!is.null(dim(theirs)) && (min(dim(theirs)) == 1L || is.null(dim(ours)))) {
    theirs <- as.numeric(theirs)
  }
  if (!is.null(dim(theirs)) && !identical(dim(as.matrix(ours)), dim(theirs))) {
    stop("Dimension mismatch in check ", name)
  }
  if (length(ours) != length(theirs)) {
    stop("Length mismatch in check ", name, ": ", length(ours), " vs ",
         length(theirs), " (", paste(signif(theirs, 4), collapse = ", "), ")")
  }
  err <- max(abs(as.numeric(ours) - as.numeric(theirs))) /
    max(1, max(abs(as.numeric(theirs))))
  rows[[length(rows) + 1L]] <<- data.frame(check = name, rel_error = err)
}
for (i in seq_along(ref$care)) {
  p <- lapply(prob$care[i, ], function(v) m(v[[1]]))
  add(paste0("care[", i, "]"), ctrl_care(p$A, p$B, p$Q, p$R)$P, ref$care[[i]])
}
for (i in seq_along(ref$dare)) {
  p <- lapply(prob$dare[i, ], function(v) m(v[[1]]))
  add(paste0("dare[", i, "]"), ctrl_dare(p$A, p$B, p$Q, p$R)$P, ref$dare[[i]])
}
lp <- lapply(prob$lyap[1, ], function(v) m(v[[1]]))
# both use A X + X A' + Q = 0
add("lyap", ctrl_lyap(lp$A, lp$Q), ref$lyap[[1]])
dp <- lapply(prob$dlyap[1, ], function(v) m(v[[1]]))
add("dlyap", ctrl_lyap(dp$A, dp$Q, discrete = TRUE), ref$dlyap[[1]])
for (i in seq_along(ref$hinfnorm)) {
  p <- lapply(prob$hinfnorm[i, ], function(v) m(v[[1]]))
  add(paste0("hinfnorm[", i, "]"),
      ctrl_norm(ctrl_ss(p$A, p$B, p$C, p$D), "inf"), ref$hinfnorm[[i]])
}
p <- lapply(prob$h2norm[1, ], function(v) m(v[[1]]))
add("h2norm", ctrl_norm(ctrl_ss(p$A, p$B, p$C, p$D), "2"), ref$h2norm[[1]])
p <- lapply(prob$hsv[1, ], function(v) m(v[[1]]))
add("hsv", ctrl_hsvd(ctrl_ss(p$A, p$B, p$C, p$D)), ref$hsv[[1]])
p <- prob$c2d[1, ]
d <- ctrl_c2d(ctrl_ss(m(p$A[[1]]), m(p$B[[1]]), m(p$C[[1]]), m(p$D[[1]])), p$dt)
add("c2d A", d$A, ref$c2d[[1]]$A)
add("c2d B", d$B, ref$c2d[[1]]$B)
p <- prob$place[1, ]
K <- ctrl_place(m(p$A[[1]]), m(p$B[[1]]), unlist(p$poles))
add("place", sort(Re(eigen(m(p$A[[1]]) - m(p$B[[1]]) %*% K)$values)), ref$place[[1]])
h <- ref$hinfsyn[[1]]
if (!is.null(h$error)) stop("The python-control hinfsyn reference failed: ", h$error)
{
  p <- prob$hinfsyn[1, ]
  P <- ctrl_ss(m(p$A[[1]]), m(p$B[[1]]), m(p$C[[1]]), m(p$D[[1]]))
  syn <- ctrl_hinfsyn(P, p$nmeas, p$ncon, tol = 1e-6)
  # both are bisection results; compare the achieved levels loosely
  add("hinfsyn gamma", syn$gamma, h$gamma)
}
res <- do.call(rbind, rows)
res$pass <- res$rel_error < ifelse(grepl("hinfsyn", res$check), 1e-3, 1e-7)
print(res, row.names = FALSE)
cat("python-control", ref$versions$control, "\n")
if (!all(res$pass)) stop("Reference comparison failed.")
