\name{print.slidevector}
\alias{print.slidevector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print Method for the Slide-vector Model
}
\description{
Print method for objects of class slidevector.
}
\usage{

\method{print}{slidevector}(x, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Object of class slidevector
}
  \item{\dots}{
Ignored
}
}

\examples{
dis <- matrix(c(1,2,3,4,5,6,2,8,9,3), nrow = 5, ncol = 2)
a <- rbind(dis,dis+1.5)
test <- as.matrix(dist(a))[1:5,6:10]
v <- slidevector(test, 2, 250, .001)
print(v)
}
