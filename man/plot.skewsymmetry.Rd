\name{plot.skewsymmetry}
\alias{plot.skewsymmetry}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting method for the skew-symmetric part of an asymmetric matrix
}
\description{
This plotting method provides a multidimensional representation of skew-symmmetry based on the singlular value decomposition (SVD). The properties of the SVD of a skew-symmetric matrix and a novel way of plotting skew-symmetry was given by Gower (1977). The singular vectors come in pairs with equal singular values. These diagrams are not interpreted by comparing distances between point, but by comparing areas formed by two points and the origin. These singular vectors span a plane, and the area of the triangle between two points and the origin represents skew-symmetry. The signs of the skew-symmetry is modelled by a direction in the plane.
}
\usage{
\method{plot}{skewsymmetry}(x, plot.plane = 1, yplus = 0, xlab, ylab, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class skewsymmetry
}
  \item{plot.plane}{
Integer indicating which plane to plot
}
  \item{yplus}{
Offset for the labels above the object points
}
  \item{xlab}{
Label for the x-axis
}
  \item{ylab}{
Label for the y-axis
}
  \item{\dots}{
Further plot arguments
}
}
\references{
Gower, J.C. (1977) The analysis of asymmetry and orthogonality. In: Recent Developments in Statistics ( J. Barra, F. Brodeau, G. Romier & B. van Cutsem, Eds.), 109-123. North Holland, Amsterdam.
}
\keyword{Methods}

