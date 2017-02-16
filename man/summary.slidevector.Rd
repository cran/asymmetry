\name{summary.slidevector}
\alias{summary.slidevector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary method for the slide-vector model
}
\description{
This function extracts the configuration matrix and the slide-vector from the object.
}
\usage{
\method{summary}{slidevector}(object, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class \code{slidevector}
}
  \item{\dots}{
Ignored
}
}

\examples{
data(Englishtowns)
v <- slidevector(Englishtowns,dim = 2,itmax = 250, eps = .001)
summary(v)

}
