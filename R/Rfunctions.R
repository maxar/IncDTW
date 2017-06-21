dtw <- function(Q, C, ws = NULL, return_diffM = FALSE){
   # require(GCM)
   # wrapper function for the C++ implementation of
   # an calculation of the global cost matrix + direction matrix
   # for dynamic time warping
   #
   # Q ... one dimensional vector, query time series
   # C ... one dimensional vector, time series to be fitted to C, and for which values are constantly observed

   #--- initial checking
   n <- length(Q)
   m <- length(C)
   if(!is.null(ws)){
      if(abs(n-m)>ws){
         stop("window size is too small, no warping path can be found")
      }
      if(ws>=m){
         stop("window size is bigger than m, no warping path can be found")
      }
   }


   #--- preparation
   mQ <- matrix(Q, ncol = m, nrow = n, byrow = F)
   mC <- matrix(C, ncol = m, nrow = n, byrow = T)
   diffM <- mQ - mC
   cm <- abs(mQ-mC)
   rm(list = c("mQ", "mC"))

   #--- calculation in C++
   if(is.null(ws)){
      ret <- GCM_cpp(cm)
   } else {
      ret <- GCM_Sakoe_cpp(cm, ws)
      # print("gcm calc done")
   }
   tmp <- BACKTRACK_cpp(ret$dm)
   ret <- c(ret, list(ii = rev(tmp$ii),
                      jj = rev(tmp$jj),
                      wp = rev(tmp$wp)))

   # ret$ii <- rev(ret$ii)#not really necessary
   # ret$jj <- rev(ret$jj)

   if(return_diffM) ret <- c(ret, list(diffM = diffM))

   return(ret)
}


idtw <- function(Q, C, newO, gcm, dm, ws = NULL, return_diffM = FALSE){
   # require(GCM)
   # wrapper function for the C++ implementation of
   # an incremental calculation of the global cost matrix + direction matrix
   # for dynamic time warping
   #
   # Q ... one dimensional vector, query time series
   # C ... one dimensional vector, time series to be fitted to C, and for which values are constantly observed
   # newO ... one dimensional vector, new observation, to be appended to C
   # gcm ... global cost matrix, output from dtw(Q,C)
   # dm ... direction matrix, output from dtw(Q,C)

   #--- initial checking
   n <- length(Q)
   m <- length(C)
   o <- length(newO)
   m2 <- o+m
   if(!is.null(ws)){
      if(abs(n-m2)>ws){
         stop("window size is too small, no warping path can be found")
      }
      if(ws>=m){
         stop("window size is bigger than m, no warping path can be found")
      }
   }

   #--- preparation
   C <- c(C, newO)
   # gcm <- cbind(gcm, matrix(0, nrow = nrow(gcm), ncol = o))
   # dm  <- cbind(dm, matrix(0, nrow = nrow(dm), ncol = o))
   gcm <- cbind(gcm, matrix(NA, nrow = n, ncol = o))
   dm  <- cbind(dm, matrix(NA, nrow = n, ncol = o))
   mQ <- matrix(Q, ncol = o, nrow = n, byrow = F)
   mC <- matrix(newO, ncol = o, nrow = n, byrow = T)
   diffM <- matrix((mQ - mC), ncol = length(newO))
   cm <- abs(diffM)

   #--- calculation in C++
   if(is.null(ws)){
      ret <- IGCM_cpp(gcmN = gcm, dmN = dm, cmN = cm)
   } else {
      ret <- IGCM_Sakoe_cpp(gcmN = gcm, dmN = dm, cmN = cm, ws = ws)
   }
   tmp <- BACKTRACK_cpp(ret$dm)
   ret <- c(ret, list(ii = rev(tmp$ii),
                      jj = rev(tmp$jj),
                      wp = rev(tmp$wp)))

   if(return_diffM) ret <- c(ret, list(diffM = diffM))
   return(ret)
}
