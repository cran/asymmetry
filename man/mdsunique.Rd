\name{mdsunique}
\alias{mdsunique}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MDS Model with Unique Dimensions
}
\description{
This MDS model has both common and unique dimensions. The common dimensions are shared by all other objects, whereas unique dimension apply to one object. A unique dimension has a non zero value for only one object, the coordinates for the other objects are zero. There are as many unique dimensions as there are objects. An asymmetric version of this model has two sets of unique dimensions: one for the rows and one for the columns. The distance in this model is defined as:
\deqn{d_{ij}(X)=\sqrt{\sum_{s=1}^p (x_{is}-x_{js})^2 + r_{i}^{2}+c_{j}^{2}}}.
}
\usage{
mdsunique(data, dim = 2, verbose = FALSE, itmax = 125, eps = 1e-12)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A data matrix
%%     ~~Describe \code{data} here~~
}
  \item{dim}{ The number of dimensions
%%     ~~Describe \code{dim} here~~
}
  \item{verbose}{ If true, prints the iteration history to screen
%%     ~~Describe \code{verbose} here~~
}
  \item{itmax}{   Maximum number of iterations
%%     ~~Describe \code{itmax} here~~
}
  \item{eps}{  Convergence criterion for Stress
%%     ~~Describe \code{eps} here~~
}
}
\value{
  \item{ndim }{The number of dimensions}
  \item{fulldim }{The number of dimensions of the full model}
  \item{stress }{The raw stress for this model}
  \item{confi}{Returns the configuration matrix of shared dimensions of this multidimensional scaling model}
  \item{X}{Returns the configuration matrix of the full model consisting of shared and unique dimensions}
  \item{niter}{The number of iterations for the algorithm to converge}
  \item{nobs}{The number of observations in this model}
  \item{resid}{A matrix with raw residuals}

}
\author{ Berrie Zielman
}

\examples{
\dontrun{
data("studentmigration")
mm<-studentmigration
mm[mm==0]<-.5          # replace zeroes by a small number
mm <- -log(mm/sum(mm)) # convert similarities to dissimilarities
v<-mdsunique(mm,dim = 2,itmax = 2100,verbose=TRUE, eps = .0000000001)
plot(v, yplus = .3)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{MDS}% use one of  RShowDoc("KEYWORDS")