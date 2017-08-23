dtw <- function(Q, C, ws = NULL, return_diffM = FALSE, return_QC = FALSE){
   # require(GCM)
   # wrapper function for the C++ implementation of
   # an calculation of the global cost matrix + direction matrix
   # for dynamic time warping
   #
   # Q ... one dimensional vector, query time series,
   #           or matrix of differences/ costs (diffM, cm)
   # C ... one dimensional vector, time series to be fitted to C, 
   #           and for which values are constantly observed.
   #           if Q is not a vector C needs to be one 
   #           of the following strings ('diffM', 'cm')
  
   if(is.character(C)){
      if(C == "diffM"){
         cm <- abs(Q)
         n <- nrow(cm)
         m <- ncol(cm)
         diffM <- Q
         rm(list=c("Q"))
         
      } else if(C == "cm"){
         cm <- Q
         n <- nrow(cm)
         m <- ncol(cm)
         diffM <- NA
      } else{
         stop("C needs to be a vector or one of the two strings: 'diffM' for difference Matrix or 'cm' for cost Matrix")
      }
   } else {
      n <- length(Q)
      m <- length(C)
   }
  
   #--- initial checking
   if(!is.null(ws)){
      if(abs(n-m)>ws){
         stop("window size is too small, no warping path can be found")
      }
   }

   #--- preparation
   if(!is.character(C)){
      mQ <- matrix(Q, ncol = m, nrow = n, byrow = F)
      mC <- matrix(C, ncol = m, nrow = n, byrow = T)
      diffM <- mQ - mC
      cm <- abs(diffM)
      rm(list = c("mQ", "mC"))
   }
   
   #--- calculation in C++
   if(is.null(ws)){
      ret <- GCM_cpp(cm)
   } else {
      ret <- GCM_Sakoe_cpp(cm, ws)
   }
   tmp <- BACKTRACK_cpp(ret$dm)
   ret <- c(ret, list(ii = rev(tmp$ii),
                      jj = rev(tmp$jj),
                      wp = rev(tmp$wp)))

   if(return_diffM) ret <- c(ret, list(diffM = diffM))
   if(return_QC) ret <- c(ret, list(Q = Q, C = C))
   class(ret) <- append(class(ret), "idtw")
   return(ret)
}


idtw <- function(Q, C, newO, gcm, dm, ws = NULL, 
                 return_diffM = FALSE, return_QC = FALSE){
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
   }

   #--- preparation
   C <- c(C, newO)
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
   if(return_QC) ret <- c(ret, list(Q = Q, C = C))
   class(ret) <- append(class(ret), "idtw")
   return(ret)
}


dec_dm <- function(dm, Ndec){
  # wrapper function for the C++ implementation of
  # an decremental calculation of the warping path subject to a given
  # direction matrix for dynamic time warping
  # 
  # dm ... direction matrix, output from dtw(Q,C)
  # Ndec ... integer, number of observations (columns) to be reduced
  
  Nnew <- ncol(dm) - Ndec
  if(Nnew == 1){
     warning("Nnew = 1, calculation is maybe meaningless")
     ret <- list(ii = nrow(dm):1,
                 jj = rep(1, nrow(dm)),
                 wp = rep(3, nrow(dm)-1) )
  } else {
     tmp <- BACKTRACK_cpp(dm[, 1:Nnew])
     ret <- list(ii = rev(tmp$ii),
                 jj = rev(tmp$jj),
                 wp = rev(tmp$wp))
  }
  return(ret)
}



plot.idtw <- function(x, type = "QC", ...) {

   pm <- pmatch(type, c("QC", "warp"))
   if(is.na(pm)) stop("The Parameter type needs to be one of c('QC', 'warp')")

   switch(pm,
          plotQC(x, ...),
          plotWarp(x, ...)
   )
}

## an alias
plot_idtw <- plot.idtw;



plotQC <- function(x, Q = NULL, C = NULL, ...){
   # Q <- sin(1:10)#+rnorm(20)
   # # C <- sin(1:30)+rnorm(30)
   # C <- sin(-2:10)#+rnorm(15)
   # # C <- c(Q, cumsum(rnorm(20))) + rnorm(120, 0, 0.5)
   # x <- dtw(Q = Q, C = C,  return_diffM = FALSE, return_QC = T)
   
   if(is.null(Q) & !is.null(x$Q)) Q <- x$Q
   if(is.null(C) & !is.null(x$C)) C <- x$C
   if(is.null(Q) | is.null(C)){
      stop("Q and C either need to be returned from dtw() or idtw(), or defined manually")
   }
   tmp1 <- data.frame(val = c(Q, C),
                      x = c(1:length(Q), 1:length(C)),
                      id = c(rep("Q", length(Q)),
                             rep("C", length(C))
                      )
   )
   
   tmp2 <- data.frame(x = c(x$ii, x$jj), 
                         y = c(Q[x$ii], C[x$jj]),
                         id = rep(1:length(x$ii), 2))
   
   gg1 <- ggplot(tmp1, ...) +
      geom_line(aes_string(x = 'x', y = 'val', group = 'id', col = 'id')) + 
      xlab("Index") + ylab("Value") +
      geom_point(aes_string(x = 'x', y = 'val', col = 'id')) +
      geom_line(data = tmp2, aes_string(x = 'x', y = 'y', group = 'id'), lty = 2) +
      guides(col = guide_legend(title = NULL))
   ret <- gg1
   
   return(ret)
}


plotWarp <- function(x, Q = NULL, C = NULL, ...){
   if(is.null(Q) & !is.null(x$Q)) Q <- x$Q
   if(is.null(C) & !is.null(x$C)) C <- x$C
   if(is.null(Q) | is.null(C)){
      stop("Q and C either need to be returned from dtw() or idtw(), or defined manually")
   }
   
   tmp3 <- data.frame(ii=x$ii, jj = x$jj)
   tmp1 <- data.frame(val = c(Q, C),
                      x = c(1:length(Q), 1:length(C)),
                      id = c(rep("Q", length(Q)),
                             rep("C", length(C))
                      )
   )
   text_yQ <- mean(Q)
   text_yC <- mean(C)
   text_xwp <- length(Q)
   text_ywp <- length(C)
   
   gg_wp <- ggplot(tmp3, ...) +
      geom_line(aes_string(x='ii', y='jj')) +
      geom_point(aes_string(x = 'ii', y = 'jj')) + 
      scale_x_continuous(position = "top", name = "", 
                         breaks = scales::pretty_breaks(n=5),
                         minor_breaks = NULL) +
      scale_y_continuous(trans = "reverse", name = "", 
                         breaks = scales::pretty_breaks(n=5),
                         minor_breaks = NULL) +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())+
         geom_text(label="Warping Path: ii", aes(y=1, x = text_xwp), hjust = 1, vjust = 1) +
         geom_text(label="Warping Path: jj", aes(x=1, y = text_ywp), hjust = 0, vjust = 1, angle = 90)

   gg_Q <- ggplot(tmp1[tmp1$id == "Q", ]) + 
      geom_line(aes_string(x = 'x', y = 'val')) + ylab("")+
      geom_point(aes_string(x = 'x', y = 'val')) +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_blank()) +
      theme(panel.grid.major.y =  element_line(colour = "#eaeaea"))+
      theme(panel.grid.minor.y =  element_line(colour = "#eaeaea"))+
      scale_x_continuous( breaks = scales::pretty_breaks(n=5),
                          minor_breaks = NULL) +
      geom_text(label="Q", aes(x=1, y=text_yQ))
   
   gg_C <- ggplot(tmp1[tmp1$id == "C", ]) + 
      geom_line(aes_string(x = 'x', y = 'val')) + 
      geom_point(aes_string(x = 'x', y = 'val')) +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
   gg_C <- gg_C + coord_flip() + 
      scale_x_continuous(position = "top", trans = "reverse", name="", 
                         breaks = scales::pretty_breaks(n=5), 
                         minor_breaks = NULL) +
      scale_y_continuous(name="", position="top")  +
      theme(panel.grid.major.x =  element_line(colour = "#eaeaea"))+
      theme(panel.grid.minor.x =  element_line(colour = "#eaeaea"))+
      geom_text(label="C", aes(x=1, y=text_yC))
   
   
   gg2 <- gridExtra::grid.arrange(gg_wp, gg_C, gg_Q,
                                  layout_matrix = rbind(c(1,1,1,2), 
                                                        c(1,1,1,2),
                                                        c(1,1,1,2),
                                                        c(3,3,3,NA)))
   ret <- gg2
   return(ret)
}

