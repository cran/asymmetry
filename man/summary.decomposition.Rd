\name{summary.decomposition}
\alias{summary.decomposition}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary method of the decomposition
}
\description{
Prints a decomposition of the sum of squares of an asymmetric matrix. The first column gives the sum of squares, whereas the second column gives the percentages of the two components. This decomposition can be applied to data, but also to the residuals of a fitted model.}
\usage{
\method{summary}{decomposition}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class \code{decomposition}
}
  \item{\dots}{
further parameters
}
}
\author{
Berrie Zielman
}


\examples{
data(Englishtowns)
q <- decomposition(Englishtowns)
summary(q)
}
