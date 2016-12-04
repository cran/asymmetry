\name{slidevector}
\alias{slidevector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The slide-vector model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The slide vector model is a scaling model for asymmetric data. A distance model is fitted to the symmetric part of the data whereas the asymmetric part of the data is represented by projections of the coordinates onto the slide-vector. The slide-vector points in the direction of large asymmetries in the data. The distance is modified in such a way that the distance between two points that are parallel  to the slide-vector is larger in the direction of this vector. The distance is smaller in the opposite direction. If the line connecting two points is perpendicular to the slide-vector the difference between the two projections is zero. In this case the distance between the two points is symmetric. The algorithm for fitting this model is derived from the majorization approach to multidimensional scaling.
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
The slide-vector model is a special case of the unfolding model. The stress is minimized by a version of the SMACOF algorithm. The main output are the configuration of points and the slide-vector.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
 \item{ndim }{the number of dimensions}
 \item{stress }{the raw stress for this model}
 \item{confi}{returns the configuration matrix of this multidimensional scaling model}
  \item{niter}{the number of iterations for the algorithm to converge}
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
