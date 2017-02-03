#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector inf_mat_match(NumericVector left, NumericVector right, NumericVector elp, int maxlength) {

  NumericVector x(maxlength, 0.0);
  int nrow = left.size();

  int upto = max(right);
  for(unsigned int i = 0; i <= upto; i++) {
    for(unsigned int j = 0; j < nrow; j++) {
      if((left[j] <= i) && (i < right[j]))
        x[i] += elp[j];
    }
  }

  return x;

}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R

*/
