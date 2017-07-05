mdsunique <-
  function(data,dim=2,verbose=FALSE,itmax=125,eps=1e-12){
    if (sum(data<0) > 0) stop("data for the slide-vector model should be positive")
    #check for number of rows equals number of columns
    if (nrow(data)!=ncol(data)) stop("the same number of rows and columns are expected")
    hjw <- .jnew("asymmetry/unfolding")
    dis <- .jcall(hjw, "[[D","unfoldingmodel",.jarray(data,dispatch=TRUE),as.integer(dim),verbose,as.integer(itmax),as.double(eps))
    ndim <- .jcall(hjw,"I","getDimension")
    nobs <- .jcall(hjw,"I","getNobs")
    niter <- .jcall(hjw,"I","getNiter")
    stress <- .jcall(hjw,"D","getStress")
#    res <- t(sapply(.jcall(hjw, "[[D","getConfiguration"),.jevalArray))
    X = t(sapply(dis,.jevalArray))
    mat <- X[1:nobs,] #  include unicity
    if(!is.null(rownames(data)))
      rownames(mat)<-rownames(data)
    colnames(mat) <- paste("D",1:(dim(mat)[2]),sep="")
    pred <- as.matrix(dist(X,upper = TRUE))
    pred <- pred[1:nobs,(nobs+1):(2*nobs)]
    resid <- data -pred
    confi <- mat[,1:dim]
    result<-list(ndim=dim,stress=stress,X=X,confi= confi,resid = resid,fulldim=ndim,
                 niter=niter,nobs=nobs)
    class(result)<-"mdsunique"
    return(result)
  }
