\name{idtw}
\alias{idtw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Incremental Dynamic Time Warping
}
\description{
Wrapper function for the incremental designed C++ calculations of the global cost matrix and the BACKTRACK_cpp that finds the cheapest warping path.
}
\usage{
idtw(Q, C, newO, gcm, dm, ws = NULL, return_diffM = FALSE, return_QC = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{Q}{one dimensional numeric vector}
\item{C}{one dimensional numeric vector}
\item{newO}{one dimensional numeric vector of new observations to be appended to C}
\item{gcm}{global cost matrix, output from dtw(Q=Q, C=C, ws=ws)}
\item{dm}{direction matrix, output from dtw(Q=Q, C=C, ws=ws)}
\item{ws}{integer, describes the window size for the sakoe chiba window. If NULL, then no window is applied. (default = NULL)}
\item{return_diffM}{boolean, if TRUE then the Matrix of differences (not the absolute value) is returned. (default = FALSE)}
\item{return_QC}{
boolean, if TRUE then the input vectors Q and C are appended to the returned list. This is useful for the \code{\link{plot.idtw}} function.  (default = FALSE)}
}
\details{
The dynamic time warping distance is the element in the last row and last column of the global cost matrix.
}
\value{
\item{gcm}{global cost matrix}
\item{dm}{direction matrix (3=up, 1=diagonal, 2=left)}
\item{wp}{warping path}
\item{ii}{indices of C of the optimal path}
\item{jj}{indices of Q of the optimal path}
\item{diffM}{Matrix of differences}
\item{Q}{input Q}
\item{C}{input C}
}
\references{
Sakoe, H.; Chiba, S., Dynamic programming algorithm optimization for spoken word recognition, Acoustics, Speech, and Signal Processing [see also IEEE Transactions on Signal Processing], IEEE Transactions on , vol.26, no.1, pp. 43-49, Feb 1978.
http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1163055
}
\author{
Maximilian Leodolter
}

\examples{
Q <- cumsum(rnorm(100))
C <- Q[11:100] + rnorm(90, 0, 0.5)
tmp <- dtw(Q=Q, C=C, ws = 15, return_diffM = FALSE) # the ordinary calculation
# new observation
newObs <-  c(2,3)
# the incremental step
idtw(Q, C, newO = newObs, gcm=tmp$gcm, dm=tmp$dm, ws = 15, return_diffM = FALSE) 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }% use one of  RShowDoc("KEYWORDS")
\keyword{ ts }% __ONLY ONE__ keyword per line
