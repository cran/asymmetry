\name{slidevector}
\alias{slidevector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The slide-vector model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function fits the slide vector model to asymmetric data.
}
\usage{
slidevector(data, dim = 2, itmax = 125, eps = 1e-12)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{A matrix with asymmetric relationship data}
\item{dim}{The number of dimensions for this model}
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
v <- slidevector(Englishtowns,2,250,.0001)
plot(v)



}
