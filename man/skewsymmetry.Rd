\name{skewsymmetry}
\alias{skewsymmetry}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Decompose an Asymmetric Matrix into Symmetric and Skew-symmetric Components
}
\description{
The decomposition of an asymmetric matrix into a symmetric matrix and a skew-symmetric matrix is an elementary result from mathematics that is the cornerstone of this package. The decomposition into a skew-symmetric and a symmetric component is written as: \eqn{Q=S+A}, where \eqn{Q} is an asymmetric matrix, \eqn{S} is a symmetric matrix, and \eqn{A} is a skew-symmetric matrix. This decomposition provides a justification for separate analyses of \eqn{S} and \eqn{A}. This decompositon is a usefull tool for data analysis, and a second application is to the study of an asymmetric matrix of residuals, obtaind after fitting a model.
}
\usage{
skewsymmetry(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An asymmetric matrix
}
}
\value{
  \item{S }{The symmetric part of the matrix}
  \item{A }{The skew-symmetric part of the matrix}
  \item{linear}{The linear part of the skew-symmetric matrix}
  \item{sv}{The singular vectors of the skew-symmetric matrix}
  \item{nobj}{The number of objects}

}
\author{
Berrie Zielman
}
\examples{
data("Englishtowns")
Q <- skewsymmetry(Englishtowns)
# the skew-symmetric part
Q$A
}
\seealso{
 \code{\link{plot.skewsymmetry}}
}

\keyword{Methods}
