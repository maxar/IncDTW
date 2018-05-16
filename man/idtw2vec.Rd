\name{idtw2vec}
\alias{idtw2vec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dynamic Time Warping
}
\description{
Wrapper function for the C++ calculations of vector based fast calculation of DTW. Allows incremental calculation for new observations and sakoe chiba warping window.
}
\usage{
idtw2vec(Q, newObs, gcm_lc = NULL, nC = NULL, ws = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Q}{one dimensional vector, query time series.}
\item{newObs}{one dimensional vector of new observations}
\item{gcm_lc}{vector, last column of global cost matrix of previous calculation. If NULL, then DTW is calculated and the last column and last row are returned to start upcoming incremental calculations. (default = NULL)}
\item{nC}{
integer, is the length of the original time series C, of which newObs are the new observations. Length of time series C exclusive new observations, such that length(c(C,newObs)) = length(newObs) + nC. Necessary if ws is not NULL.  (default = NULL)
}
\item{ws}{
integer, describes the window size for the sakoe chiba window. If NULL, then no window is applied. (default = NULL)
}
}


\details{
no matrices are allocated, no matrices are returned
}
\value{
\item{distance}{the DTW distance}
\item{gcm_lc_new}{the last column of the new global cost matrix}
\item{gcm_lr}{the last row of the new global cost matrix}
}
\references{
Sakoe, H.; Chiba, S., Dynamic programming algorithm optimization for spoken word recognition, Acoustics, Speech, and Signal Processing [see also IEEE Transactions on Signal Processing], IEEE Transactions on , vol.26, no.1, pp. 43-49, Feb 1978. http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1163055
}
\author{
Maximilian Leodolter
}

\examples{
Q <- cumsum(rnorm(100))
C <- Q[11:100] + rnorm(90, 0, 0.5)

# initial calculation
res0 <- idtw2vec(Q=Q, newObs = C, gcm_lc = NULL)

# incremental calculation for new observations
res1 <- idtw2vec(Q, newObs = rnorm(10), gcm_lc = res0$gcm_lc_new)


}

\keyword{ classif }% use one of  RShowDoc("KEYWORDS")
\keyword{ cluster }% use one of  RShowDoc("KEYWORDS")
\keyword{ ts }% __ONLY ONE__ keyword per line

