\name{dtw2vec}
\alias{dtw2vec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dynamic Time Warping
}
\description{
Wrapper function for the C++ calculations of vector based fast calculation of DTW. Allows early abandoning and sakoe chiba warping window.
}
\usage{
dtw2vec(Q, C, ws = NULL, threshold = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Q}{one dimensional vector, query time series.}
\item{C}{one dimensional vector, query time series.}
\item{ws}{
integer, describes the window size for the sakoe chiba window. If NULL, then no window is applied. (default = NULL)
}
\item{threshold}{numeric, the threshold for early abandoning. In the calculation of the global cost matrix a possible path stops as soon as the threshold is reached. Facilitates faster calculations incase of low threshold. (default = FALSE)
}
}

\details{
no matrices are allocated, no matrices are returned
}
\value{
\item{distance}{the DTW distance}
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
dtw2vec(Q = Q, C = C)
dtw2vec(Q = Q, C = C, ws = 30)
dtw2vec(Q = Q, C = C, threshold = 100)
dtw2vec(Q = Q, C = C, ws = 30, threshold = 100)

}
\keyword{ classif }% use one of  RShowDoc("KEYWORDS")
\keyword{ cluster }% use one of  RShowDoc("KEYWORDS")
\keyword{ ts }% __ONLY ONE__ keyword per line

