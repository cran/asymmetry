\name{asymscal}
\alias{asymscal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
asymscal
}
\description{
This function fits the asymscal model proposed by Young (1975). Asymmetry in the data arises by weighting the dimensions of a multidimensional scaling configuration. When a subject compares object \emph{i} to \emph{j} he or she may use different weights when comparing object \emph{j} to \emph{i}.
}
\usage{
asymscal(data, ndim = 2, itmax = 10000, eps = 1e-10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The data
}
  \item{ndim}{
The number of dimensions
}
  \item{itmax}{
The maximum number of iterations
}
  \item{eps}{
Convergence criterion for Stress
}
}
\details{
This function exploits a connection between the INDSCAL model and the asymscal model. This method inherits the methods for plotting an printing from the \code{smacofIndDiff} in the smacof package. Basically, the asymscal takes two steps. First, this function sets up the appropriate dissimilarity and missing data structure for a three-way multidimensional scaling model, then a call to the method \code{smacofIndDiff} in the imported package smacof is made. After correcting for the normalization applied to the data by \code{smacofIndDiff}, the results can be displayed and plotted by the methods in the package smacof.
}
\value{
  \item{delta}{Observed dissimilarities}
  \item{obsdiss}{List of observed dissimilarities, normalized}
  \item{confdiss}{List of configuration dissimilarities}
  \item{conf}{List of matrices of final configurations}
  \item{gspace}{Joint configurations aka group stimulus space}
  \item{cweights}{Configuration weights}
  \item{stress}{Stress-1 value}
  \item{resmat}{Matrix with squared residuals}
  \item{rss}{Residual sum-of-squares}
  \item{spp}{Stress per point}
  \item{ndim}{Number of dimensions}
  \item{model}{Type of smacof model}
  \item{niter}{Number of iterations}
  \item{nobj}{Number of objects}
}
\references{
Young, F. W. (1975). An asymmetric Euclidean model for multi-process asymmetric data. Paper presented at the U.S.-Japan Seminar on Multidimensional scaling, San Diego, U.S.A.
}
\author{
Berrie Zielman
}
\note{
The original algorithm for fitting the asymscal model fits squared distances. This function fits distances and not squared distances. The configuration matrix is normalized to one.
}

\section{Warning }{Although my experience with fitting this model in the majorization framework is very limited, it appears to converge slowly. I recommend keeping the convergence criteria in the example}

\examples{
\dontrun{
data("asymscalexample")
t<-asymscal(data,ndim=2,itmax=10000,eps=1e-10)
t$cweights
round(t$cweights,3)
plot(t, plot.type = "confplot")
plot(t, plot.type = "bubbleplot")
plot(t, plot.type = "stressplot")
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")