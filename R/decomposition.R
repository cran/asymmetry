decomposition <- function(X)
{
  S <-.5*(t(X)+X)
  A <-.5*(X - t(X))
  lk <- list(S = S,A = A)
  class(lk) <- "decomposition"
  lk
}

