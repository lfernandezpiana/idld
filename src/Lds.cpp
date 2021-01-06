// LambdaBeta c++ version

#include <Rcpp.h>
#include <cmath> 
#include <iostream>

using namespace Rcpp;

double cppAbs(double x) {
  if (x>0) {
    return x;
  } else {
    return (-1)*x;
  }
}

double cppEmpirica(double z, NumericVector data) {
  int n = data.size();
  int e = 0;
  for (int i=0; i<n; i++) { 
    if ( z >= data[i] ) {
      e++;
    }
    else {
      e = e;
    }
  }  
  double emp = e;
  return emp/n;
}


double cppLambda (double z,NumericVector data, float beta)
{
  int n = data.size();
  int k = n*beta;
  double w;
  NumericVector orden = clone(data);
  for (int i=0; i<n; i++) {
  w = z - data[i];
  orden[i] = cppAbs(w);
  }
  std::sort(orden.begin(), orden.end());
  return orden[k-1];
}

// [[Rcpp::export]]
double cppLds(double z,NumericVector data,float beta) {
  int n = data.size();
  double l = cppLambda(z,data,beta);
  double sup = cppEmpirica(z+l,data) - cppEmpirica(z,data);
  double inf = cppEmpirica(z,data) - cppEmpirica(z-l,data);
  double lds = (2/pow(beta,2))*sup*inf; 
  return lds;
}

// [[Rcpp::export]]
double cppLdaux(NumericVector W,NumericMatrix proyecciones,float beta) {
  int m = W.size();
  double depth = 0;
  double depthAux; 
  for (int i=0; i < m; i++) {
    depthAux = cppLds(W[i],proyecciones.column(i),beta);
    depth = depthAux + depth;
  }
  double depthAux2 = depth/m;
  return depthAux2;
}

// [[Rcpp::export]]
NumericVector cppLdaux2(NumericMatrix W, NumericMatrix proyecciones, float beta, int q) {
  NumericVector lidd (q); 
  for (int i=0; i<q; i++){
  lidd[i] = cppLdaux(W.column(i), proyecciones, beta);   
  }
  return lidd;
}
  
