print.mdsunique<-
  function(x,...)
  {
    cat("Dimensions:              ")
    cat(x$ndim)
    cat("\n")
    cat("Number of objects:       ")
    cat(x$nobj)
    cat("\n")
    cat("Number of iterations:    ")
    cat(x$niter)
    cat("\n")
    cat("Stress:                   ")
    cat(x$stress)
    cat("\n")
  }
