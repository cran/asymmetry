slidevector <-
function(data,dim=2,verbose=FALSE,itmax=125,eps=1e-12){
  if (sum(data<0) > 0) stop("data for the slide-vector model should be positive")
  if (nrow(data)!=ncol(data)) stop("the same number of rows and columns are expected")
  hjw <- .jnew("asymmetry/slidevector")     # create instance of slidevector class
  dis <- .jcall(hjw, "[[D","slidevectormodel",.jarray(data,dispatch=TRUE),as.integer(dim),verbose,as.integer(itmax),as.double(eps))
  ndim <- .jcall(hjw,"I","getDimension")
  nobs <- .jcall(hjw,"I","getNobs")
  niter <- .jcall(hjw,"I","getNiter")
  stress <- .jcall(hjw,"D","getStress")
  X = t(sapply(dis,.jevalArray))
  mat <- X[1:nobs,]
  Z <- X[nobs+1,]
  if(!is.null(rownames(data)))
    rownames(mat)<-rownames(data)
  colnames(mat) <- paste("D",1:(dim(mat)[2]),sep="")
  a <- as.matrix(rbind(mat+matrix(Z,nrow=nobs,ncol=ndim,byrow=TRUE),mat))
  pred <- as.matrix(dist(a,upper = TRUE))
  pred <- pred[1:nobs,(nobs+1):(2*nobs)]
  resid <- data-pred
  result<-list(ndim=ndim,stress=stress,confi=mat,slvec=Z,resid = resid, niter=niter,nobs=nobs,pred=pred)
  class(result)<-"slidevector"
  return(result)
}
