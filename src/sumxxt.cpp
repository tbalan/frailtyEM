#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector sumxxt(List x, int L) {
  List xlist(x);
  int n = xlist.size();

  std::vector<double> res((L * (L+1)) / 2);


  for(int l = 0; l < n; l++) {
    // Rcout<<"at line"<<l<<std::endl;

    SEXP ll = xlist[l];
    Rcpp::NumericVector xi(ll);

    for(int j = 0; j < L; j++)
      for(int i = j; i < L; i++) {
        // Rcout<<"adding "<<xi[i]<<"*"<<xi[j]<<"="<<xi[i] * xi[j]<<" to "<<(i * (i+1))/2 + j<<" so ("<<i<<", "<<j<<")"<<std::endl;
        res[(i * (i+1))/2 + j] += xi[i] * xi[j];
      }
  }

  return wrap(res);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
# This function is used to calculate if you have a list of vectors of the same length
# and you want the sum of all the x%*%x transposed
# n <- 1
# L <- 3
# set.seed(1)
# x <- matrix(rnorm(n * L), nrow = n, ncol = L)
# xx <- apply(x, 1, list)
# xxx <- lapply(xx, function(x) x[[1]])
#
# xxx
#
# m <- matrix(0, L, L)
# a <- sumxxt(xxx, L)
# m[upper.tri(m, diag = TRUE)] <- a
#
# m3 <- Reduce("+", lapply(xxx, function(vec) vec %*% t(vec)))
#
# m2 <- m + t(m) - diag(diag(m))
# all.equal(m2, m3)
# all.equal(m[upper.tri(m)] , m3[upper.tri(m3)])


*/
