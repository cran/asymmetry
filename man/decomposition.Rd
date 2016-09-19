\name{decomposition}
\alias{decomposition}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Decompose an asymmetric matrix
}
\description{
The decomposition of an asymmetric matrix into a symmetric matrix and a skew-symmetric matrix is an elementary result from mathematics that is the cornerstone of this package. The decomposition into a skew-symmetric and a symmetric component is written as: \eqn{Q=S+A}, where \eqn{Q} is an asymmetric matrix, \eqn{S} is a symmetric matrix, and \eqn{A} is a skew-symmetric matrix.

}
\usage{
decomposition(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
An asymmetric matrix
}
}
\value{
  \item{S }{The symmetric part of the matrix}
  \item{A }{The skew-symmetric part of the matrix}

}
\author{
Berrie Zielman
}
\examples{
data("Englishtowns")
Q <- decomposition(Englishtowns)
# the skew-symmetric part
Q$A

}
