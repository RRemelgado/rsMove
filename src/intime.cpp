# include <Rcpp.h>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]

NumericMatrix intime(NumericMatrix x, NumericVector ij, NumericVector oj, NumericVector z) {

  int nrow = x.nrow();
  int ncol = x.ncol();
  int nrec = oj.size();

  NumericMatrix res(nrow,nrec);

  for (int r=0; r < nrow; r++) {

    for (int i=0; i < nrec; i++) {

      int si = ncol+1;
      for (int j=0; j < ncol; j++) {if((ij[j] == oj[i]) & (LogicalVector::is_na(x(r,j)) == 0)) {si = j;}}

      if (si > ncol) {

        // find indices (before)
        int bi = ncol+1;
        int mn = oj[i]-z[0];
        int mx = oj[i];
        for (int j=0; j < ncol; j++) {if((LogicalVector::is_na(x(r,j)) == 0) & (ij[j] >= mn) & (ij[j] < mx)) {bi = j;}}

        // find indices (after)
        int ai = ncol+1;
        mn = oj[i];
        mx = oj[i]+z[1];
        for (int j=0; j < ncol; j++) {if((LogicalVector::is_na(x(r,j)) == 0) & (ai > ncol) & (ij[j] > mn) & (ij[j] < mx)) {ai = j;}}

        if ((bi < ncol) & (ai < ncol)) {

          res(r,i) = x(r,bi) + (oj[i]-ij[ai]) * ((x(r,ai)-x(r,bi)) / (ij[ai]-ij[bi]));

        } else {res(r,i) = NumericVector::get_na();}
      } else {res(r,i) = x(r,si);}

    }

  }

  return(res);

}
