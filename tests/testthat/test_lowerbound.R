context("manual testing")

test_that("lowerbound smaller than dtw distance", {
   # test it with the function cpp_get_tube
   # which is not exported to an R fucntion but also 
   # available after devtools::load_all(".")
   # => the final counter should be 0!
  skip("")
   
   
   
   ws <- 10
   nh <- 100
   nx <- 100
   h <- cumsum(rnorm(nh))
   x <- cumsum(rnorm(nx))
   tube <- cpp_get_tube(h, ws);
   j0 <- 0
   jsup <- 0+nh
   get_lb(tube, x, 0, jsup)
   dtw2vec(h, x[(j0+1):(jsup)], ws=ws, dist_method = "norm1", step_pattern = "symmetric1")$distance
   
   
   
   plot(h, type="l")
   lines(tube[,1], col="blue")
   lines(tube[,2], col="red")
   
   
   
   
   counter <- 0
   for(i in 1:1000){
      myseed <- sample(10^4, size = 1)
      set.seed(myseed)
      
      nh <- 100
      nx <- 100
      ws <- sample(nh, size = 1)
      h <- cumsum(rnorm(nh))
      x <- cumsum(rnorm(nx))
      tube <- cpp_get_tube(h, ws);
      j0 <- 0
      jsup <- 0+nh
      lb <- get_lb(tube, x, 0, jsup)
      d <- dtw2vec(h, x[(j0+1):(jsup)], ws=ws, dist_method = "norm1", step_pattern = "symmetric1")$distance
      
      if(lb > d){
         print(myseed)
         counter <- counter +1
      }
   }
   print(counter)
   
   
})

