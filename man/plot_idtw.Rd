\name{plot_idtw}
\alias{plot_idtw}
\alias{plot.idtw}
\alias{plotQC}
\alias{plotWarp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the results from Dynamic Time Warping
}
\description{
plot function for the output of \code{\link{dtw}} and \code{\link{idtw}}}
\usage{
\method{plot}{idtw}(x, type = "QC", ...)
# an alias for plot_idtw
plot_idtw(x, type="QC", ...)

plotQC(x, Q, C, ...)

plotWarp(x, Q, C, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{output from dtw(Q, C)}
 \item{Q}{one dimensional numeric vector}
\item{C}{one dimensional numeric vector}
\item{type}{character, one of c('QC', 'warp')}
\item{...}{Other arguments passed on to methods. Not currently used.}
}
\details{
The plot function visualizes the time warp and the alignment of the two time series.}
\author{
Maximilian Leodolter
}

\examples{
Q <- cumsum(rnorm(100))
C <- Q[11:100] + rnorm(90, 0, 0.5)
tmp <- dtw(Q = Q, C = C, ws = 15, return_wp = TRUE, return_QC = TRUE) # the ordinary calculation
plotQC(tmp, Q, C)
plotWarp(tmp, Q, C)
plot(tmp, type = 'warp')
plot(tmp, Q, C, type = 'QC')
plot(tmp, type = 'QC')# since return_QC = TRUE, this is also possible
}
