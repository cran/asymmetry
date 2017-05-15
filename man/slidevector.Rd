\name{slidevector}
\alias{slidevector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The slide-vector model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The slide-vector model is a multidimensional scaling model for asymmetric data. Here, a distance model is fitted to the symmetric part of the data whereas the asymmetric part of the data is represented by projections of the coordinates onto the slide-vector. The slide-vector points in the direction of large asymmetries in the data. The interpretation of asymmetry in this model is aided by the use of projections of points onto the slide-vector. The distance from i to j is larger if the point i has a higher projection onto the slide-vector than the distance from j to i. If the line connecting two points is perpendicular to the slide-vector the difference between the two projections is zero. In this case the distance between the two points is symmetric. The algorithm for fitting this model is derived from the majorization approach to multidimensional scaling.
}
\usage{
slidevector(data, dim = 2, verbose = FALSE, itmax = 125, eps = 1e-12)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{An asymmetric matrix}
\item{dim}{The number of dimensions for this model}
\item{verbose}{Print the history of iterations}
\item{itmax}{The maximum number of iterations}
\item{eps}{The convergence criterion for the algorithm}
}
\details{
The slide-vector model is a special case of the unfolding model. Therefore, the algorithm for fitting this model is a constrained unfolding model. The coordinates of the objects are calculated by minimizing a least squares loss function. This loss function is called stress in the multidimensional scaling literature. The stress is minimized by a version of the SMACOF algorithm. The main output are the configuration of points and the slide-vector.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
 \item{ndim }{The number of dimensions}
 \item{stress }{The raw stress for this model}
 \item{confi}{Returns the configuration matrix of this multidimensional scaling model}
  \item{niter}{The number of iterations for the algorithm to converge}
  \item{nobs}{the number of observations}
  \item{resid}{A matrix with raw residuals}
  \item{slvec }{Coordinates of the slide-vector}

}
\references{
Zielman, B., and Heiser, W. J. (1993), The analysis of asymmetry by a slide-vector, Psychometrika, 58, 101-114.
}
\author{
Berrie Zielman
}

\seealso{
 \code{\link{plot.slidevector}}
}
\examples{
## asymmetric distances between English towns
data(Englishtowns)
v <- slidevector(Englishtowns,dim = 2,itmax = 250, eps = .001)
plot(v)
}
\keyword{models}
