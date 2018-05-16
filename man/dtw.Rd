\name{dtw}
\alias{dtw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dynamic Time Warping
}
\description{
Wrapper function for the C++ calculations of the global cost matrix and the BACKTRACK_cpp that finds the cheapest warping path.
}
\usage{
dtw(Q, C, ws = NULL, return_cm = FALSE,
                     return_diffM = FALSE,
                     return_wp = FALSE,
                     return_diffp = FALSE,
                     return_QC = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Q}{one dimensional vector, query time series. Alternatively Q can also be a matrix of differences/ costs (diffM, cm).
}
\item{C}{
one dimensional vector, time series to be fitted to C, and for which values are constantly observed. If Q is not a vector but a matrix, then C needs to be one of the following strings ('diffM', 'cm').
}
\item{ws}{
integer, describes the window size for the sakoe chiba window. If NULL, then no window is applied. (default = NULL)
}
\item{return_cm}{boolean, if TRUE then the Matrix of costs (the absolute value) is returned. (default = FALSE)
}
\item{return_diffM}{
boolean, if TRUE then the Matrix of differences (not the absolute value) is returned. (default = FALSE)
}
\item{return_wp}{
boolean, if TRUE then the warping path is returned. (default = FALSE) If return_diffp == TRUE, then return_wp is set to TRUE as well.
}
\item{return_diffp}{
boolean, if TRUE then the path of differences (not the absolute value) is returned. (default = FALSE)
}
\item{return_QC}{
boolean, if TRUE then the input vectors Q and C are appended to the returned list. This is useful for the \code{\link{plot.idtw}} function.  (default = FALSE)
}
}
\details{
The dynamic time warping distance is the element in the last row and last column of the global cost matrix.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{distance}{the DTW distance, that is the element of the last row and last column of gcm}
\item{gcm}{global cost matrix}
\item{dm}{direction matrix (3=up, 1=diagonal, 2=left)}
\item{wp}{warping path}
\item{ii}{indices of C of the optimal path}
\item{jj}{indices of Q of the optimal path}
\item{cm}{Matrix of costs}
\item{diffM}{Matrix of differences}
\item{diffp}{path of differences}
\item{Q}{input Q}
\item{C}{input C}
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
dtw(Q = Q, C = C, ws = 15, return_diffM = FALSE)

# compare different input variations
dtw_base <- dtw(Q = Q, C = C, ws = 15, return_diffM = TRUE)
dtw_diffM <- dtw(Q = dtw_base$diffM, C = "diffM", ws = 15, return_diffM = TRUE)
dtw_cm <- dtw(Q = abs(dtw_base$diffM), C = "cm", ws = 15, return_diffM = TRUE)

identical(dtw_base$gcm, dtw_cm$gcm)
identical(dtw_base$gcm, dtw_diffM$gcm)

# of course no diffM is returned in the 'cm'-case
dtw_cm$diffM
}
\keyword{ cluster }% use one of  RShowDoc("KEYWORDS")
\keyword{ ts }% __ONLY ONE__ keyword per line

