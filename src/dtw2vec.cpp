#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
//using namespace stdlib; 

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
// [[Rcpp::plugins(cpp11)]]


double mymin(double x, double y)
{
   // z > nan for z != nan is required by C the standard
   int xnan = isnan(x), ynan = isnan(y);
   if(xnan || ynan) {
      if(xnan && !ynan) return y;
      if(!xnan && ynan) return x;
      return x;
   }
   return std::min(x,y);
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
double cpp_dtw2vec (NumericVector x, NumericVector y)
{
   
   int nx = x.size();
   int ny = y.size();
   double * p1 = new double [nx];
   double * p2 = new double [nx];
   double * ptmp; 
   double ret;
   
   // first column
   *p1 = abs(x[0]-y[0]);
   for(int i=1; i<nx; i++){
      p1[i] = abs(x[i]-y[0]) + p1[i-1];
   }
   

   for(int j=1; j < ny; j++){
      *p2 = abs(x[0]-y[j]) + *(p1);
      
      for(int i=1; i<nx; i++){
         *(p2+i) = abs(x[i]-y[j]) + mymin(*(p2+i-1), mymin(*(p1+i), *(p1+i-1)));
      }
      ptmp=p1;
      p1 = p2;
      p2 = ptmp;
   }
   
   ret = *(p1+nx-1);//p1[nx-1]
   delete[] p1;
   delete[] p2;
   
   return (ret);
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
List cpp_dtw2vec_inc (NumericVector x, NumericVector newObs, NumericVector gcm_lc)
{
   // x ... time series that is fixed
   // y ... time series of new observations, ONLY new observations to be appended
   // gcm_lc ... last column of old GCM
   // TODO...
   int nx = x.size();
   int nnewObs = newObs.size();
   double * p1 = new double [nx];
   double * p2 = new double [nx];
   double * ptmp; 
   double mynan;
   NumericVector gcm_lr(nnewObs);
   NumericVector gcm_lc_new(nx);
   mynan = std::numeric_limits<double>::quiet_NaN();
   
   if(nx != gcm_lc.size()){
      return mynan;
   }
   
   // first column
   for(int i=0; i<nx; i++){
      p1[i] = gcm_lc[i];
   }
   
   for(int j=0; j < nnewObs; j++){
      *p2 = abs(x[0]-newObs[j]) + *(p1);
      
      for(int i=1; i<nx; i++){
         *(p2+i) = abs(x[i]-newObs[j]) + mymin(*(p2+i-1), mymin(*(p1+i), *(p1+i-1)));
      }
      gcm_lr[j] = *(p2+nx-1);
      ptmp=p1;
      p1 = p2;
      p2 = ptmp;
   }
   
   for(int i=0; i<nx; i++){
      gcm_lc_new[i] = *(p1+i);
   }
   
   
   List ret;
   ret["gcm_lr"] = gcm_lr;
   ret["gcm_lc_new"] = gcm_lc_new;
   ret["distance"] = *(p1+nx-1);
   delete[] p1;
   delete[] p2;
   return ret ;
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
List cpp_dtw2vec_inc_ws (NumericVector x, NumericVector newObs, NumericVector gcm_lc,
                         int ws, int ny)
{
   // x ... time series that is fixed
   // newObs ... time series of new observations, ONLY new observations to be appended
   // gcm_lc ... last column of old GCM
   // int ws ... window size
   // int ny ... length of time series y exclusive new observations, 
   //             length(c(y,newObs)) = length(newObs) + ny
   
   int nx = x.size();
   int nnewObs = newObs.size();
   int iBegin = 0;
   int iEnd = 0;
   double * p1 = new double [nx];
   double * p2 = new double [nx];
   double * ptmp; 
   double mynan;
   NumericVector gcm_lr(nnewObs);
   NumericVector gcm_lc_new(nx);
   mynan = std::numeric_limits<double>::quiet_NaN();
   
   if(nx != gcm_lc.size()){
      return mynan;
   }
   
   // first column
   for(int i=0; i<nx; i++){
      p1[i] = gcm_lc[i];
      p2[i] = mynan;//initialize b with NAN
   }
   
   for(int j=0; j < nnewObs; j++){
      iBegin = ny+j-ws;
      if(iBegin <= 0){
         *p2 = abs(x[0]-newObs[j]) + *(p1);
         iBegin = 1;   
      }else if (iBegin == 1){
         *(p2+iBegin -1) = mynan;//must not be available for finding the cheapest path
      }else{
         *(p2+iBegin -1) = mynan;//must not be available for finding the cheapest path
         *(p2+iBegin -2) = mynan;//must not be available for finding the cheapest path
      }
      
      iEnd   = ny+j+ws+1;
      if(iEnd >= nx){
         iEnd = nx;   
      }else{
         *(p1+iEnd) = mynan;//must not be available for finding the cheapest path
      }
      
      for (int i = iBegin; i < iEnd; i++){
         *(p2+i) = abs(x[i]-newObs[j]) + mymin(*(p2+i-1), mymin(*(p1+i), *(p1+i-1)));
      }
      gcm_lr[j] = *(p2+nx-1);
      ptmp=p1;
      p1 = p2;
      p2 = ptmp;
   }
   
   for(int i=0; i<nx; i++){
      gcm_lc_new[i] = *(p1+i);
   }
   
   List ret;
   ret["gcm_lr"] = gcm_lr;
   ret["gcm_lc_new"] = gcm_lc_new;
   ret["distance"] = *(p1+nx-1);
   delete[] p1;
   delete[] p2;
   
   return ret ;
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
double cpp_dtw2vec_ws (NumericVector x, NumericVector y, int ws)
{
   int iBegin = 0;
   int iEnd = 0;
   
   int nx = x.size();
   int ny = y.size();
   double * p1 = new double [nx];
   double * p2 = new double [nx];
   double * ptmp; 
   double mynan;
   double ret;
   mynan = std::numeric_limits<double>::quiet_NaN();
   
   
   //initialize a and b with NAN
   for(int i=0; i < nx; i++){
      p1[i] = mynan;
      p2[i] = mynan;
   }
   
   // first column
   p1[0] = abs(x[0]-y[0]);
   iEnd   = std::min(nx, ws);
   for(int i=1; i < iEnd; i++){
      p1[i] = abs(x[i]-y[0]) + p1[i-1];
   }
   
   
   for(int j=1; j < ny; j++){
      iBegin = j-ws;
      if(iBegin <= 0){
         *p2 = abs(x[0]-y[j]) + *(p1);
         iBegin = 1;   
      }else if (iBegin == 1){
         *(p2+iBegin -1) = mynan;//must not be available for finding the cheapest path
      }else{
         *(p2+iBegin -1) = mynan;//must not be available for finding the cheapest path
         *(p2+iBegin -2) = mynan;//must not be available for finding the cheapest path
      }
      
      
      iEnd   = j+ws+1;
      if(iEnd >= nx){
         iEnd = nx;   
      }else{
         *(p1+iEnd) = mynan;//must not be available for finding the cheapest path
      }
      
      for (int i = iBegin; i < iEnd; i++){
         *(p2+i) = abs(x[i]-y[j]) + mymin(*(p2+i-1), mymin(*(p1+i), *(p1+i-1)));
      }
      ptmp=p1;
      p1 = p2;
      p2 = ptmp;
   }
   ret = *(p1+nx-1);//p1[nx-1]
   delete[] p1;
   delete[] p2;
   
   return (ret);
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
double cpp_dtw2vec_ea (NumericVector x, NumericVector y, double threshold)
{
   
   int nx = x.size();
   int ny = y.size();
   double * p1 = new double [nx];
   double * p2 = new double [nx];
   double * ptmp; 
   double mynan;
   double ret;
   int nanCounter;
   mynan = std::numeric_limits<double>::quiet_NaN();
   
   // first column
   p1[0] = abs(x[0]-y[0]);
   if(p1[0] > threshold)  return(mynan); 
   
   for(int i=1; i<nx; i++){
      p1[i] = abs(x[i]-y[0]) + p1[i-1];
      if(p1[i] > threshold) p1[i] = mynan;
   }

   for(int j=1; j < ny; j++){
      
      nanCounter = 0;
      *p2 = abs(x[0]-y[j]) + *(p1);
      if(*(p2) > threshold){
         *(p2) = mynan;
         nanCounter ++;
      }
      
      for(int i=1; i<nx; i++){
         *(p2+i) = abs(x[i]-y[j]) + mymin(*(p2+i-1), mymin(*(p1+i), *(p1+i-1)));
         if((*(p2+i) > threshold) | (*(p2+i) != *(p2+i))){
            *(p2+i) = mynan;
            nanCounter ++;
         } 
      }
      if(nanCounter == nx) return(mynan);
      ptmp=p1;
      p1 = p2;
      p2 = ptmp;
   }
   
   ret = *(p1+nx-1);//p1[nx-1]
   delete[] p1;
   delete[] p2;
   
   return (ret);
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
double cpp_dtw2vec_ws_ea (NumericVector x, NumericVector y, int ws,double threshold)
{
   int iBegin = 0;
   int iEnd = 0;
   
   int nx = x.size();
   int ny = y.size();
   double * p1 = new double [nx];
   double * p2 = new double [nx];
   double * ptmp; 
   double mynan;
   double ret;
   int nanCounter;
   mynan = std::numeric_limits<double>::quiet_NaN();
   
   
   //initialize a and b with NAN
   for(int i=0; i < nx; i++){
      p1[i] = mynan;
      p2[i] = mynan;
   }
   
   // first column
   p1[0] = abs(x[0]-y[0]);
   if(p1[0] > threshold)  return(mynan); 
   
   iEnd   = std::min(nx, ws);
   for(int i=1; i < iEnd; i++){
      p1[i] = abs(x[i]-y[0]) + p1[i-1];
      if(p1[i] > threshold) p1[i] = mynan;
   }
   
  
   
   for(int j=1; j < ny; j++){
      nanCounter = 0;
      iBegin = j-ws;
      if(iBegin <= 0){
         *p2 = abs(x[0]-y[j]) + *(p1);
         if(*(p2) > threshold){
            *(p2) = mynan;
            nanCounter ++;
         }
         iBegin = 1;   
      }else if (iBegin == 1){
         nanCounter = 1;
         *(p2+iBegin -1) = mynan;//must not be available for finding the cheapest path
      }else{
         nanCounter = iBegin;
         *(p2+iBegin -1) = mynan;//must not be available for finding the cheapest path
         *(p2+iBegin -2) = mynan;//must not be available for finding the cheapest path
      }
      
      iEnd   = j+ws+1;
      if(iEnd >= nx){
         iEnd = nx;   
      }else{
         *(p1+iEnd) = mynan;//must not be available for finding the cheapest path
      }
      
      for (int i = iBegin; i < iEnd; i++){
         *(p2+i) = abs(x[i]-y[j]) + mymin(*(p2+i-1), mymin(*(p1+i), *(p1+i-1)));
         if((*(p2+i) > threshold) | (*(p2+i) != *(p2+i))){
            *(p2+i) = mynan;
            nanCounter ++;
         }
      }
      
      if(nanCounter == nx) return(mynan);
      ptmp=p1;
      p1 = p2;
      p2 = ptmp;
   }
   ret = *(p1+nx-1);//p1[nx-1]
   delete[] p1;
   delete[] p2;
   
   return (ret);
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


// [[Rcpp::export]]
double cpp_dtw2vec_v32 (NumericVector x, NumericVector y)
{
   
   int nx = x.size();
   int ny = y.size();
   int i;
   double * p1 = new double [nx];
   double * p2 = new double [nx];
   double * ptmp; 
   double ret;
   
   // first column
   p1[0] = abs(x[0]-y[0]);
   for(i=1; i<nx; i++){
      p1[i] = abs(x[i]-y[0]) + p1[i-1];
   }
   
   for(int j=1; j < ny; j++){
      *p2 = abs(x[0]-y[j]) + *(p1);
      
      for( i=1; i < nx-32; i+= 32){
         *(p2+i)   = abs(x[i]-y[j])   + mymin(*(p2+i-1), mymin(*(p1+i)  , *(p1+i-1)));
         *(p2+i+1) = abs(x[i+1]-y[j]) + mymin(*(p2+i)  , mymin(*(p1+i+1), *(p1+i)));
         *(p2+i+2) = abs(x[i+2]-y[j]) + mymin(*(p2+i+1), mymin(*(p1+i+2), *(p1+i+1)));
         *(p2+i+3) = abs(x[i+3]-y[j]) + mymin(*(p2+i+2), mymin(*(p1+i+3), *(p1+i+2)));
         *(p2+i+4) = abs(x[i+4]-y[j]) + mymin(*(p2+i+3), mymin(*(p1+i+4), *(p1+i+3)));
         *(p2+i+5) = abs(x[i+5]-y[j]) + mymin(*(p2+i+4), mymin(*(p1+i+5), *(p1+i+4)));
         *(p2+i+6) = abs(x[i+6]-y[j]) + mymin(*(p2+i+5), mymin(*(p1+i+6), *(p1+i+5)));
         *(p2+i+7) = abs(x[i+7]-y[j]) + mymin(*(p2+i+6), mymin(*(p1+i+7), *(p1+i+6)));
         *(p2+i+8) = abs(x[i+8]-y[j]) + mymin(*(p2+i+7), mymin(*(p1+i+8), *(p1+i+7)));
         *(p2+i+9) = abs(x[i+9]-y[j]) + mymin(*(p2+i+8), mymin(*(p1+i+9), *(p1+i+8)));
         *(p2+i+10) = abs(x[i+10]-y[j]) + mymin(*(p2+i+9), mymin(*(p1+i+10), *(p1+i+9)));
         *(p2+i+11) = abs(x[i+11]-y[j]) + mymin(*(p2+i+10), mymin(*(p1+i+11), *(p1+i+10)));
         *(p2+i+12) = abs(x[i+12]-y[j]) + mymin(*(p2+i+11), mymin(*(p1+i+12), *(p1+i+11)));
         *(p2+i+13) = abs(x[i+13]-y[j]) + mymin(*(p2+i+12), mymin(*(p1+i+13), *(p1+i+12)));
         *(p2+i+14) = abs(x[i+14]-y[j]) + mymin(*(p2+i+13), mymin(*(p1+i+14), *(p1+i+13)));
         *(p2+i+15) = abs(x[i+15]-y[j]) + mymin(*(p2+i+14), mymin(*(p1+i+15), *(p1+i+14)));
         
         *(p2+i+16) = abs(x[i+16]-y[j]) + mymin(*(p2+i+15), mymin(*(p1+i+16), *(p1+i+15)));
         *(p2+i+17) = abs(x[i+17]-y[j]) + mymin(*(p2+i+16), mymin(*(p1+i+17), *(p1+i+16)));
         *(p2+i+18) = abs(x[i+18]-y[j]) + mymin(*(p2+i+17), mymin(*(p1+i+18), *(p1+i+17)));
         *(p2+i+19) = abs(x[i+19]-y[j]) + mymin(*(p2+i+18), mymin(*(p1+i+19), *(p1+i+18)));
         *(p2+i+20) = abs(x[i+20]-y[j]) + mymin(*(p2+i+19), mymin(*(p1+i+20), *(p1+i+19)));
         *(p2+i+21) = abs(x[i+21]-y[j]) + mymin(*(p2+i+20), mymin(*(p1+i+21), *(p1+i+20)));
         *(p2+i+22) = abs(x[i+22]-y[j]) + mymin(*(p2+i+21), mymin(*(p1+i+22), *(p1+i+21)));
         *(p2+i+23) = abs(x[i+23]-y[j]) + mymin(*(p2+i+22), mymin(*(p1+i+23), *(p1+i+22)));
         
         *(p2+i+24) = abs(x[i+24]-y[j]) + mymin(*(p2+i+23), mymin(*(p1+i+24), *(p1+i+23)));
         *(p2+i+25) = abs(x[i+25]-y[j]) + mymin(*(p2+i+24), mymin(*(p1+i+25), *(p1+i+24)));
         *(p2+i+26) = abs(x[i+26]-y[j]) + mymin(*(p2+i+25), mymin(*(p1+i+26), *(p1+i+25)));
         *(p2+i+27) = abs(x[i+27]-y[j]) + mymin(*(p2+i+26), mymin(*(p1+i+27), *(p1+i+26)));
         *(p2+i+28) = abs(x[i+28]-y[j]) + mymin(*(p2+i+27), mymin(*(p1+i+28), *(p1+i+27)));
         *(p2+i+29) = abs(x[i+29]-y[j]) + mymin(*(p2+i+28), mymin(*(p1+i+29), *(p1+i+28)));
         *(p2+i+30) = abs(x[i+30]-y[j]) + mymin(*(p2+i+29), mymin(*(p1+i+30), *(p1+i+29)));
         *(p2+i+31) = abs(x[i+31]-y[j]) + mymin(*(p2+i+30), mymin(*(p1+i+31), *(p1+i+30)));
         
      }
      for (; i < nx; i++) {
         *(p2+i)   = abs(x[i]-y[j])   + mymin(*(p2+i-1), mymin(*(p1+i)  , *(p1+i-1)));
      }
      ptmp=p1;
      p1 = p2;
      p2 = ptmp;
   }
   ret = *(p1+nx-1);//p1[nx-1]
   delete[] p1;
   delete[] p2;
   
   return (ret);
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


