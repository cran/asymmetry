summary.mdsunique <-
function(object,...)
{
  cat("\n")
  cat("Configurations:\n")
  xda <- cbind(object$confi,object$row,object$col)
  nam <- colnames(object$confi)
  names <- c(nam,"Row","Col")
  colnames(xda) <- names
  print(round(xda,4))
}
