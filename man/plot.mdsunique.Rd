\name{plot.mdsunique}
\alias{plot.mdsunique}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Configuration Plot for mdsunique
}
\description{
Method for a two-dimensional plot of the model. Available rownames are plotted as labels above the points.
}
\usage{
\method{plot}{mdsunique}(x, plot.dim = c(1, 2), yplus = 0, xlab, ylab, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Object of class \code{mdsunique}
}
  \item{plot.dim}{
A vector with dimensions to be plotted
}
  \item{yplus}{
  Parameter to adjust the vertical position of the label
  }
  \item{xlab}{
Label of x-axis.
}
  \item{ylab}{
Label of y-axis.
}
  \item{\dots}{Further plot arguments: see \code{\link[graphics]{plot}} for detailed information.}
}
\examples{
\dontrun{
data("studentmigration")
mm<-studentmigration
mm[mm==0]<-.5          # replace zeroes by a small number
mm <- -log(mm/sum(mm)) # convert similarities to dissimilarities
v<-mdsunique(mm,ndim = 2,itmax = 2100,verbose=TRUE, eps = .0000000001)
plot(v, yplus = .3)
}
}
\keyword{Methods}
