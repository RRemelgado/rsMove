#include <Rcpp.h>
using namespace Rcpp;

//' @title runmean2
//'
//' @description NA-sensitive running mean.
//' @param x A \emph{Numeric} vector.
//' @param y A \emph{numeric} element.
//' @return A vector.
//' @details {Applies a running mean over \emph{x} with a window size defined by \emph{y}.
//' Missing values are ignored and/or filled during the compution.}
//' @export
// [[Rcpp::export]]

NumericVector runmean2(NumericVector x, int y) {

  // get dimensions of x
  int ne = x.size();

  // apply running mean filter
  NumericVector res = x;

  for (int j=1;j<ne-1; j++) {

    // index of values used to estimate the mean
    NumericVector index = NumericVector::create(j-y,j,j+y);
    NumericVector v = x[index];
    res[j] = mean(na_omit(v));

  }

  return(res);

}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
runmean2(c(0.1,0.2,0.3,NA,0.5,0.4,NA,0.2), 1)
  */
