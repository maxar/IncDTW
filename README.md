# IncDTW - Incremental Dynamic Time Warping


### What is this github repo for?
This is a developement version of the package available on CRAN, so some features and functions are available here before I commit to CRAN.

### Get Stared
```R
#install.packages("devtools")
library(devtools)
install_github("maxar/IncDTW")

library(IncDTW)
help("IncDTW-package")
help("dtw2vec")
```

### Package Description
The Dynamic Time Warping (DTW) distance measure for time series allows non-linear alignments of 
time series to match similar patterns in time series of different lengths and or different speeds. 

IncDTW is characterized by 
(1) the incremental calculation of DTW (reduces runtime complexity to a linear level for 
updating the DTW distance) - especially for life data streams or subsequence matching, 
(2) the vector based implementation of DTW which is faster because 
no matrices are allocated (reduces the space complexity from a quadratic to a linear level in the number of observations) - for all
runtime intensive DTW computations, 
(3) the subsequence matching algorithm runDTW, that efficiently finds the k-NN to a query pattern in a long time series, and 
(4) C++ in the heart. 

For details about DTW a good start to read is the wikipedia page <https://en.wikipedia.org/wiki/Dynamic_time_warping>, where
you will find lots of references, as the original paper "Dynamic programming algorithm optimization for spoken word recognition"
by Sakoe and Chiba (1978) <https://ieeexplore.ieee.org/document/1163055>.


