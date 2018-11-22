# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export()]]
NumericMatrix intime(NumericVector x, NumericMatrix y, NumericVector z) {

  int nrow = y.nrow();
  int ncol = y.ncol();

  NumericMatrix res(nrow,ncol);

  for (int r=0; r < nrow; r++) {

    res(r,0) = y(r,0); //  index 1st element (can't interpolate)
    res(r,ncol-1) = y(r,ncol-1); //  index last element (can't interpolate)

    for (int i=1; i< ncol; i++) {

      if (LogicalVector::is_na(y(r,i))) {

        // find indices (before)
        int bi = nrow+1;
        int mx = x[i]-z[1];
        for (int j=i; j > -1; --j) {if((LogicalVector::is_na(y(r,j)) == 0) & (bi > nrow) & (x[j] > mx)) {bi = j;}}

        // find indices (after)
        int ai = nrow+1;
        mx = x[i]+z[2];
        for (int j=i+1; j < nrow; j++) {if((LogicalVector::is_na(y(r,j)) == 0) & (ai > nrow) & (x[j] > mx)) {ai = j;}}

        if ((bi < nrow) & (ai < nrow)) {

          double a = 0;
          double b = 0;

          double xsum = x[bi] + x[ai];
          double ysum = y(r,bi) + y(r,ai);
          double x2sum = pow(x[bi],2) + pow(x[ai],2);
          double xsum2 = pow(x[ai]+x[bi],2);
          double xysum = (x[ai*y(r,ai)]) + (x[ai]*x[ai]);

          a = ((ysum*x2sum) - (xsum*xysum)) / (2*x2sum-xsum2);
          b = ((xysum) - (xsum*ysum)) / (2 * x2sum - xsum2);

          res(r,i) = a + b * x[i];

        } else {res(r,i) = NumericVector::get_na();}
      } else {res(r,i) = y(r,i);}

    }

  }

  return(res);

}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
m <- matrix(0,15,15) # test matrix
m[] <- sample(1:length(m), length(m)) # assign random numbers
m[sample(1:length(m), 30)] <- NA # assign random NA's
om <- intime(c(1:ncol(m)), m, c(2,2)) # interpolate NA's
  */
