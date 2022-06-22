slidevector <-
  function(data, weight = NULL, ndim = 2, verbose = FALSE, itmax = 125, eps = 1e-12 )
  {
    if (sum(data < 0) > 0) stop("data for the slide-vector model should be non-negative")
    if (nrow(data) != ncol(data)) stop("the same number of rows and columns are expected in the data matrix")
    nobj <- nrow(data)
    if(is.null(weight))  weight <- 0*data+1
    if (sum(weight<0) > 0) stop("weight for the slide-vector model should be positive")
    if (nrow(weight) != ncol(weight)) stop("the same number of rows and columns are expected in the weight matrix")
    we <-  rbind(cbind(weight * 0, weight),cbind(t(weight), weight * 0))
    v <- diag(rowSums(we)) - we
    da <-  rbind(cbind(data*0, data), cbind(t(data), data*0))
    ztil <- matrix(runif(nobj*2*ndim),2*nobj,ndim)
    dist <- as.matrix(dist(ztil))
    ind <- dist<.000000001
    dist[ind]<-1
    BMAT <- we*da/dist
    BMAT[ind]<-0
    BMAT <- diag(rowSums(BMAT)) - BMAT #VX
    z <- rbind(diag(nobj), diag(nobj))
    z <- cbind(z,c(rep(1,nobj), rep(0,nobj)))
    zvz <-  t(z)%*%v%*%z
    nsp <- c(rep(1,nobj),0)%*%t(c(rep(1,nobj),0))/nobj
    zvz_inv <- solve(zvz+nsp)-nsp
    dist <- as.matrix(dist(ztil))
    stressold <- sum(weight*(data-dist[1:nobj,(nobj+1):(2*nobj)])^2)
    for (i in 1:itmax){
      ind <- dist<.000000001
      dist[ind]<-1
      BMAT <- we*da/dist
      BMAT[ind]<-0
      BMAT <- diag(rowSums(BMAT))-BMAT #VX
      zup <- BMAT%*%ztil
      ztil <- z%*%zvz_inv%*%t(z)%*%zup
      dist <- as.matrix(dist(ztil))
      stress <- sum(weight*(data-dist[1:nobj,(nobj+1):(2*nobj)])^2)
      if(verbose)
        print(stress)
      if(stressold - stress < eps) break
      stressold <- stress
    }
    niter <- i
    X <- zvz_inv%*%t(z)%*%zup #moet v tussen
    mat <- X[1:nobj,]
    Z <- X[nobj+1,]
    if(!is.null(rownames(data)))
      rownames(mat)<-rownames(data)
    colnames(mat) <- paste("D",1:(dim(mat)[2]),sep="")
    pred <- dist[1:nobj,(nobj+1):(2*nobj)]
    resid <- data-pred
    result<-list(ndim=ndim,stress=stress,confi = mat, slvec = Z,resid = resid, niter=niter,
                 nobs=nobj,nobj=nobj,pred=pred,model="Slide-vector model",ztil=ztil)
    class(result)<-"slidevector"
    return(result)
  }
print.slidevector <-
  function(x, ...)
  {
    cat("Dimensions:              ")
    cat(x$ndim)
    cat("\n")
    cat("Number of objects:       ")
    cat(x$nobs)
    cat("\n")
    cat("Number of iterations:    ")
    cat(x$niter)
    cat("\n")
    cat("Stress:                   ")
    cat(x$stress)
    cat("\n")
  }
summary.slidevector <-
  function(object,...)
  {
    cat("\n")
    cat("Configurations:\n")
    print(round(object$confi,4))
    cat("\n")
    cat("Slide-vector:\n")
    cat(round(object$slvec,4))
  }
