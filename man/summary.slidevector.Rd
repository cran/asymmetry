\name{summary.slidevector}
\alias{summary.slidevector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
summary method for the slide-vector model
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
An object of class \code{object} slidevector
}
  \item{\dots}{
Ignored
}
}

\examples{
data(Englishtowns)
v<-slidevector(Englishtowns,2,250,.0001)
summary(v)

}
