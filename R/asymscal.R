#
#    use weights in the next version
#

asymscal <- function(data, ndim = 2, itmax = 10000, eps = 1e-10)
{
  nrow = nrow(data)
  l <- vector("list", nrow)
  dawe <- vector("list", nrow)
  dimscal <- vector("list", nrow)
  for (i in 1:nrow)
  {
    temp <- data*0
    temp[,i] <- data[i,]
    temp[i,] <- data[i,]
    l[[i]] <- as.dist(temp)
    temp[,i] <- rep(1,nrow)
    temp[i,] <-  rep(1,nrow)
    temp[i,i] <- 0
    dawe[[i]] <- as.dist(temp)

  }
  #
  #  tweak the smacof normalization
  #
  normalizer <- lapply(l,function(x) sqrt(sum(x^2)))

  for (i in 1:nrow)
  {
    l[[i]] <- sqrt(.5*nrow*(nrow-1))*l[[i]]/normalizer[[i]]
  }
  #   end tweak function
  asy<- smacofIndDiff(l, weightmat = dawe, ndim = ndim, itmax = itmax, eps = eps, constraint = "indscal")
  wtemp<-matrix(0,nrow,ndim)
  for (i in 1:nrow)
  {
    dimscal[[i]] <- asy$cweights[[i]]*normalizer[[i]]/sqrt(.5*nrow*(nrow-1))
    asy$conf[[i]] <- asy$conf[[i]]*normalizer[[i]]/sqrt(.5*nrow*(nrow-1))
    for (j in 1:ndim)
    {
      wtemp[i,j] <- dimscal[[i]][j,j]
    }
  }
  #
  # all quantities should be recalculated to account for the normalization (confi)
  #
  sq <- apply(asy$gspace, 2, function(x) sqrt(sum(x^2)))
  wtemp <- sweep(wtemp,2,sq,"*")
  gspace2 <- sweep(asy$gspace,2,sq,"/")
  results <- list(delta = asy$delta, obsdiss = asy$obsdiss, confdiss = asy$confdiss, conf = asy$conf,
                  gspace = gspace2, cweights = wtemp, stress = asy$stress,
                  resmat = asy$resmat, rss = asy$rss, spp = asy$spp,  ndim = asy$ndim,
                  model = "ASYMSCAL", niter = asy$niter, nobj = asy$nobj, call = match.call())
  class(results) <- "smacofID"
  results
}

