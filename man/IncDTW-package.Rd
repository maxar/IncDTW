\name{IncDTW-package}
\alias{IncDTW-package}
%\alias{idtw}
\docType{package}
\title{
  Incremental Dynamic Time Warping
}
\description{
  Implements incremental calculation of the DTW ( Dynamic Time Warping) distance of two vectors, which is specifically useful for life data streams. Further the calculation of the global cost matrix is implemented in C++ to be faster. The Sakoe Chiba band is also implemented. The calculation of DTW is less functional then the one of dtw(), however much faster.
}
\details{
  Any feedback is highly appreciated!
}
\author{
Maximilian Leodolter

Maintainer: Maximilian Leodolter <maximilian.leodolter@gmail.com>
}
\references{
  Dynamic programming algorithm optimization for spoken word recognition by Sakoe and Chiba published in 1978 (DOI:10.1109/TASSP.1978.1163055)
}
\keyword{ DTW }
\seealso{
  https://github.com/maxar/IncDTW
}
