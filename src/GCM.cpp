#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List GCM_Sakoe_cpp(Rcpp::NumericMatrix cM, int ws){
   // Global Cost Matrix
   int n = cM.nrow();
   int m = cM.ncol();
   int iBegin = 0;
   int iEnd = 0;
   double cost = 0;

   NumericMatrix gcm(n, m);
   NumericMatrix dm(n, m);
   std::fill( gcm.begin(), gcm.end(), NumericVector::get_na() );
   std::fill( dm.begin(), dm.end(), NumericVector::get_na() );
   
   gcm(0,0) = cM(0,0);
   for(int i =1; i <= ws; i++){
      dm(i,0)=3;
      gcm(i,0) = cM(i, 0) + gcm(i - 1, 0);
   }
   for(int j =1; j <= ws; j++){
      dm(0,j)=2;
      gcm(0, j) = cM(0, j) + gcm(0, j - 1);
   }
   
   for (int j = 1; j < m; j++){
      iBegin = std::max(1, j-ws);
      iEnd   = std::min(n, j+ws+1);

      for (int i = iBegin; i < iEnd; i++){
         cost = cM(i,j);
         if(!std::isnan(gcm(i-1, j)) && !std::isnan(gcm(i, j-1))){
            
            if(gcm(i-1, j-1) <= gcm(i-1, j) && gcm(i-1, j-1) <= gcm(i, j-1)){
               gcm(i,j) = cost + gcm(i-1, j-1);
               dm(i,j) = 1;
            } else if(gcm(i-1, j) <= gcm(i-1, j-1) && gcm(i-1, j) <= gcm(i, j-1)){
               gcm(i,j) = cost + gcm(i-1, j);
               dm(i,j) = 3;
            }else{
               gcm(i,j) = cost + gcm(i, j-1);
               dm(i,j) = 2;
            }
            
         } else if(std::isnan(gcm(i-1, j)) && std::isnan(gcm(i, j-1))){
            
            gcm(i,j) = cost + gcm(i-1, j-1);
            dm(i,j) = 1;
            
         } else if (std::isnan(gcm(i-1, j))){
            
            if(gcm(i-1, j-1) <= gcm(i , j-1)){
               gcm(i,j) = cost + gcm(i-1, j-1);
               dm(i,j) = 1;
            } else{
               gcm(i,j) = cost + gcm(i, j-1);
               dm(i,j) = 2;
            }
            
         } else{// if (std::isnan(gcm(i, j-1))){
            
            if(gcm(i-1, j-1) <= gcm(i-1 , j)){
               gcm(i,j) = cost + gcm(i-1, j-1);
               dm(i,j) = 1;
            } else{
               gcm(i,j) = cost + gcm(i-1, j);
               dm(i,j) = 3;
            }
            
         }
      }
   }

   List ret;
   ret["gcm"] = gcm;
   ret["dm"] = dm;
   return ret ;
}

//################################################################################
//################################################################################


#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List IGCM_Sakoe_cpp(Rcpp::NumericMatrix gcmN, //global costmatrix with new empty columns
                    Rcpp::NumericMatrix dmN, //direction matrix with new empty columns
                    Rcpp::NumericMatrix cmN,//local cost matrix of new observations and old constant vector
                    int ws){ 
   
   
   // C has one to many new columns that are not covered by the global cost matrix gcm so far
   // indexC gives the first index of the new observations
   // 
   int n = gcmN.nrow();
   int m = gcmN.ncol();
   int j;
   int Nnew = cmN.ncol();//number of new observations
   double cost;
   int iBegin;
   int iEnd;
   
   // first row
   if((m-Nnew) < ws){
      j = m-Nnew;
      do{
         gcmN(0,j) = gcmN(0,j-1) + cmN(0,(j-m+Nnew)); 
         dmN(0,j)  = 2;   
         j = j + 1;
      } while (j<ws);
   }
   
   // remaining
   // 
   
   for (int j = (m-Nnew); j < m; j++){
      iBegin = std::max(1, j-ws);
      iEnd   = std::min(n, j+ws+1);
      
      for (int i = iBegin; i < iEnd; i++){
         
         
         // for (int j = (m-Nnew); j < m; j++){
         //    for (int i = 1; i < n; i++){
         cost = cmN(i,(j-m+Nnew));
         
         if(!std::isnan(gcmN(i-1, j)) && !std::isnan(gcmN(i, j-1))){
            
            if(gcmN(i-1, j-1) <= gcmN(i-1, j) && gcmN(i-1, j-1) <= gcmN(i, j-1)){
               gcmN(i,j) = cost + gcmN(i-1, j-1);
               dmN(i,j) = 1;
            } else if(gcmN(i-1, j) <= gcmN(i-1, j-1) && gcmN(i-1, j) <= gcmN(i, j-1)){
               gcmN(i,j) = cost + gcmN(i-1, j);
               dmN(i,j) = 3;
            }else{
               gcmN(i,j) = cost + gcmN(i, j-1);
               dmN(i,j) = 2;
            }
            
         } else if(std::isnan(gcmN(i-1, j)) && std::isnan(gcmN(i, j-1))){
            
            gcmN(i,j) = cost + gcmN(i-1, j-1);
            dmN(i,j) = 1;
            
         } else if (std::isnan(gcmN(i-1, j))){
            
            if(gcmN(i-1, j-1) <= gcmN(i , j-1)){
               gcmN(i,j) = cost + gcmN(i-1, j-1);
               dmN(i,j) = 1;
            } else{
               gcmN(i,j) = cost + gcmN(i, j-1);
               dmN(i,j) = 2;
            }
            
         } else{// if (std::isnan(gcmN(i, j-1))){
            
            if(gcmN(i-1, j-1) <= gcmN(i-1 , j)){
               gcmN(i,j) = cost + gcmN(i-1, j-1);
               dmN(i,j) = 1;
            } else{
               gcmN(i,j) = cost + gcmN(i-1, j);
               dmN(i,j) = 3;
            }
            
         }
         
         /*
         if(gcmN(i-1, j-1) <= gcmN(i-1, j) && gcmN(i-1, j-1) <= gcmN(i  , j-1)){
            gcmN(i,j) = cost + gcmN(i-1, j-1);
            dmN(i,j) = 1;
         } else if(gcmN(i-1, j) <= gcmN(i-1, j-1) && gcmN(i-1, j) <= gcmN(i, j-1)){
            gcmN(i,j) = cost + gcmN(i-1, j);
            dmN(i,j) = 3;
         }else{
            gcmN(i,j) = cost + gcmN(i, j-1);
            dmN(i,j) = 2;
         }
          */
      }
   }
   
   List ret;
   ret["gcm"] = gcmN;
   ret["dm"] = dmN;
   return ret ;
}


//################################################################################
//################################################################################


#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List GCM_cpp(Rcpp::NumericMatrix cM){
   // Global Cost Matrix
   int n = cM.nrow();
   int m = cM.ncol();
   
   double cost;
   //std::vector<double> lc;//local costs
   NumericMatrix gcm(n, m);
   NumericMatrix dm(n, m);
   
   gcm(0,0) = cM(0,0);
   for(int i =1; i < n; i++){
      dm(i,0)=3;
      gcm(i,0) = cM(i, 0) + gcm(i - 1, 0);
   }
   for(int j =1; j < m; j++){
      dm(0,j)=2;
      gcm(0, j) = cM(0, j) + gcm(0, j - 1);
   }
   dm(0,0) = NAN;
   
   for (int i = 1; i < n; i++){
      for (int j = 1; j < m; j++){
         cost = cM(i,j);
         if(gcm(i-1, j-1) <= gcm(i-1, j) && gcm(i-1, j-1) <= gcm(i  , j-1)){
            gcm(i,j) = cost + gcm(i-1, j-1);
            dm(i,j) = 1;
         } else if(gcm(i-1, j) <= gcm(i-1, j-1) && gcm(i-1, j) <= gcm(i, j-1)){
            gcm(i,j) = cost + gcm(i-1, j);
            dm(i,j) = 3;
         }else{
            gcm(i,j) = cost + gcm(i, j-1);
            dm(i,j) = 2;
         }
         
      }
   }
   //List z = List::create( DTW ) ;
   List ret;
   ret["gcm"] = gcm;
   ret["dm"] = dm;
   return ret ;
}


//################################################################################
//################################################################################



#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List IGCM_cpp(Rcpp::NumericMatrix gcmN, //global costmatrix with new empty columns
          Rcpp::NumericMatrix dmN, //direction matrix with new empty columns
          Rcpp::NumericMatrix cmN){ //local cost matrix of new observations and old constant vector
          
          
   // C has one to many new columns that are not covered by the global cost matrix gcm so far
   // indexC gives the first index of the new observations
   // 
   int n = gcmN.nrow();
   int m = gcmN.ncol();
   int Nnew = cmN.ncol();//number of new observations
   double cost;
   
   for (int j = (m-Nnew); j < m; j++){
      gcmN(0,j) = gcmN(0,j-1) + cmN(0,(j-m+Nnew)); 
      dmN(0,j)  = 2;
   }
   //cost = n*m;
   
   for (int j = (m-Nnew); j < m; j++){
      for (int i = 1; i < n; i++){
         cost = cmN(i,(j-m+Nnew));
         if(gcmN(i-1, j-1) <= gcmN(i-1, j) && gcmN(i-1, j-1) <= gcmN(i  , j-1)){
            gcmN(i,j) = cost + gcmN(i-1, j-1);
            dmN(i,j) = 1;
         } else if(gcmN(i-1, j) <= gcmN(i-1, j-1) && gcmN(i-1, j) <= gcmN(i, j-1)){
            gcmN(i,j) = cost + gcmN(i-1, j);
            dmN(i,j) = 3;
         }else{
            gcmN(i,j) = cost + gcmN(i, j-1);
            dmN(i,j) = 2;
         }
      }
   }
   
   List ret;
   ret["gcm"] = gcmN;
   ret["dm"] = dmN;
   return ret ;
}


//################################################################################
//################################################################################


#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
// [[Rcpp::export]]
List BACKTRACK_cpp(Rcpp::NumericMatrix dm){//direction matrix with new empty columns
              
   int n = dm.nrow();
   int m = dm.ncol();
   int i = n;
   int j = m;
   int step;
   vector<int> ii;
   vector<int> jj;
   vector<int> wp;
   
  ii.push_back(i);
  jj.push_back(j);
   
   do{
      step = dm(i-1,j-1);
      if(step == 1){
         i = i - 1;
         j = j - 1;
      } else if ( step == 2){
         j = j - 1;
      } else if ( step == 3){
         i = i - 1;
      } else{
         i = 99;
         j = 99;
      }
      ii.push_back(i);
      jj.push_back(j);
      wp.push_back(step);
   } while (i > 1 || j > 1);
   
   // ii.push_back(1);
   // jj.push_back(1);
   
   List ret;
   ret["ii"] = ii;
   ret["wp"] = wp;
   ret["jj"] = jj;
   return ret ;
}