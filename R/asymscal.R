asymscal <- function(data, ndim = 2,  start = NULL, verbose = FALSE,
                     itmax = 125, eps = 1e-10)
  {
    if (sum(data < 0) > 0)
      stop("data for the asymscal model should be positive")
    if (nrow(data) != ncol(data))
      stop("the same number of rows and columns are expected")
    nobj <- nrow(data)
    l <- vector("list", nobj)
    dawe <- vector("list", nobj)
    dimscal <- vector("list", nobj)
    for (i in 1:nobj) {
      temp <- data * 0
      temp[, i] <- data[i,]
      temp[i,] <- data[i,]
      l[[i]] <- as.dist(temp)
      temp[, i] <- rep(1, nobj) #unit weights
      temp[i,] <- rep(1, nobj)
      temp[i, i] <- 0
      dawe[[i]] <- as.dist(temp)
    }
    if (!is.matrix(start))
      start <- "torgerson"
    #
    # tweak the smacof normalization
    #
    normalizer <- lapply(l, function(x)
      sqrt(sum(x ^ 2)))
    asy <-
      smacofIndDiff(
        l,
        weightmat = dawe,
        init = start,
        ndim = ndim,
        itmax = itmax,
        eps = eps,
        constraint = "indscal"
      )
    wtemp <- matrix(0, nobj, ndim)
    for (i in 1:nobj)
    {
      dimscal[[i]] <-
        asy$cweights[[i]] * normalizer[[i]] / sqrt(.5 * nobj * (nobj - 1))
      for (j in 1:ndim)
      {
        wtemp[i, j] <- dimscal[[i]][j, j]
      }
    }
    #
    # all quantities should be recalculated to account for the normalization (confi)
    #
    sq <- apply(asy$gspace, 2, function(x)
      sqrt(sum(x ^ 2)))
    wtemp <- sweep(wtemp, 2, sq, "*")
    gspace2 <- sweep(asy$gspace, 2, sq, "/")
    res <- data
    for (i in 1:nobj)
    {
      res[i, ] <- res[i, ] - as.matrix(dist(gspace2 %*% diag(wtemp[i, ])))[i, ]
    }
    stress <- sum(res ^ 2)
    spp <- colSums(res ^ 2) / stress
    resmat <- asy$resmat
    nobj <- asy$nobj
    resmat[diag(nobj) == 1] <- 0
    if (verbose)
      print(stress)

    results <- list(
      delta = data,
      obsdiss = asy$delta,
      gspace = gspace2,
      cweights = wtemp,
      stress = stress,
      resmat = resmat,
      rss = res,
      spp = spp,
      ndim = asy$ndim,
      model = "ASYMSCAL",
      niter = asy$niter,
      nobj = asy$nobj,
      call = match.call()
    )
    class(results) <- "smacofID"
    results
  }
