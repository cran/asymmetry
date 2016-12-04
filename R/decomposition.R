decomposition <- function(x)
{
  if (nrow(x)!=ncol(x)) stop("the same number of rows and columns are expected")
  S <-.5*(t(x)+x)
  A <-.5*(x - t(x))
  lin = colMeans(A)
  lk <- list(S = S, A = A, linear = lin)
  class(lk) <- "decomposition"
  lk
}

