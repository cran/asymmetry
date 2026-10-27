hmap <- function(x, dominance = FALSE, order_rows = TRUE,  col = cm.colors(n = 9), ...)
{
  if(is.matrix(x)) x <- .5*(x-t(x))
  if(is(x,"skewsymmetry")) x <-x$A
  if (nrow(x)!=ncol(x)) stop("the same number of rows and columns are expected")
  bin <- x*0
  bin[x>0] <- 1
  bin[x<0] <- -1
  x_bin <- bin
  rsbin<-rowSums(bin)
  if(order_rows) {
    x <- x[order(rsbin),order(rsbin)]
    x_bin <- bin[order(rsbin),order(rsbin)]
  }
  #  keuze tussen binair en original
  if(dominance)
  {
    heatmap(x_bin, symm = TRUE, Rowv = NA, col = col)
  }
  else {
    heatmap(x, symm = TRUE, Rowv = NA, col = col)  }
}
