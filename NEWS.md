### Version 1.1.2.03 (in dev)

##### Bug Fixes

- if parameter ws = Inf then R breaks. A new check function checks if ws = Inf and if TRUE, then ws is set to NULL, which is equivalent to the meaning of Inf.
- scale() returned NaN if the standrad deviation was 0. Now there is a check if the sd is smaller than the parameter threshold (default = 1e-5). If the sd is smaller than the threshold, then no scaling is performed, only shifting. Analogous it's implemented for min-max scaling. Also for rundtw().

##### New Features

- New Vignettes
- replaced the function norm() by scale(), same functionality. norm() still works, but prints a warning to be deprecated.
- replaced the arguments 'normalize' by 'scale' for the function rundtw(). See details of the function documentation.

##### other changes

- Changed the name of the data set "Walk" to "walk"
- For clarification I replaced 'norm' by 'scale' in the context of z-scaling and min-max-scaling (z-normalization and min-max-normalization). From now on the terminology should be clearer seperated from normalizing the DTW distance for the length of the time series.

### Version 1.1.2

##### Bug Fixes

- fixed the case of calculating the cost matrix with cm() for univariate time series with a self defined distance function.
- export the 'insert' functions for simulate_timewarp()


##### New Features

- new branch of wrapper functions around the new S3 class 'planedtw': initialize_plane(), increment(), decrement(), reverse(), refresh(). This set of function should make it easier to navigate in the plane of possible fits, to increase the usability of the functions idtw2vec() dtw_partial() that are called behind the scenes. Also plot() and print() methods for the class 'planedtw'.
- an improved 'lot-mode' for rundtw() -- where the parameter 'C' is a list of time series -- helps to keep the allocated storage low
- new parameter '...' for the function cm() allows to pass further arguments
- new S3 class 'rundtw' for results of the function rundtw()
- print and summary methods for the S3 classes 'idtw', 'dba' and 'rundtw'
- plot method for the S3 class 'rundtw'
- new parameter 'return_QC' for rundtw() for easier plotting
- is.class() for all S3 classes in the package

##### other changes

- revised all examples in the help files
- renamed DBA() to dba(). DBA() is still available, but deprecated. A warning is printed.


### Version 1.1.1

##### Bug Fixes

- fixed issue with initial best-sofar-value-in-window for rundtw(), too many unnecessary computations were completed within the first nh observations

##### New Features

- running z-normalization for the function rundtw(). Up to now only min-max-normalization
was implemented. The parameter 'normalize' now has 3 possible values, the former TRUE and FALSE are still possible to pass. They will be translated internally to '01' and 'none' and a warning message is printed, saying that TRUE and FALSE is deprecated.
- lot-mode: ('list-of-timeseries'-mode) the parameter 'C' for the function rundtw() can now also be a list of time series. So rundtw() can search for the kNN of a query pattern Q in a list of time series of varying lengths.
- new entry in the result vector 'counter' of the function rundtw(): counter["completed"]

##### other changes

- add-on in the description of the data sets. That it's z-normalized.


### Version 1.1.0

##### Bug Fixes

- add labels to result of dtw_dismat() and dtw_disvec()
- corrected assignment of ii and jj to Q and C in the description files

##### New Features

- new function: rundtw()
- new function: find_peaks()
- new feature: for simulate_timewarp(), the parameter preserve_length
- vector based (also incremental) implementation for existing cost matrix
- plot functions for DBA for multivariate time series

##### other changes

- change name of Vignette, the name visible online
- new Vignette, which is an extensive discussion of DTW, incremental DTW
      and sub sequence matching




### Version 1.0.5

##### Bug Fixes
- normalized dtw in dtw() for multivar time series 

##### New Features

- simulate_timewarp(): new function
- dtw2vec_cm() and idtw2vec_cm() : new functions included in dtw2vec() and idtw2vec()
