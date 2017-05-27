#include <Rcpp.h>
using namespace Rcpp;

/* this function takes a left position, right position, and elp, meant to be applied for each cluster
 * it returns a vector which, for all unique time points, contains the sum of elp gathered from the cluster.
 * i.e. say we have (0, 3) and (2,3) from 5 event time points (up to 4)
 * Then this will give a vector with
 * position 0 elp1
 * position 1 elp1
 * position 2 elp1 + elp2
 * position 3 elp1 + elp2
 * position 4 0.0
 */

// [[Rcpp::export]]
NumericVector inf_mat_match(NumericVector left, NumericVector right, NumericVector elp, int maxlength) {

  NumericVector x(maxlength, 0.0);
  unsigned int nrow = left.size();

  unsigned int upto = maxlength; // this could be the length of the out vector

  for(unsigned int i = 1; i <= upto; i++) {
    for(unsigned int j = 0; j < nrow; j++) {
      if((left[j] < i) && (i <= right[j]))
        x[i-1] += elp[j];
    }
  }

  return x;

}

// this one is faster I think but something is still not right with the margins of the looping
// //[[Rcpp::export]]
// NumericVector inf_mat_match(NumericVector left, NumericVector right, NumericVector summand, int ntimes) {
//
//   NumericVector x(ntimes, 0.0);
//
//   unsigned int nrow_max = left.size();
//   for(unsigned int nrow = 0; nrow < nrow_max; nrow++) {
//     // Rcout<<"at row..."<<nrow<<std::endl;
//     // Rcout<<"edges: left "<<left[nrow]<<" and right "<<right[nrow]<<std::endl;
//     for(int tp = left[nrow]; tp <= right[nrow]; tp++) {
//       // Rcout<<"x["<<tp<<"] += summand["<<nrow<<"]"<<std::endl;
//       x[tp] += summand[nrow];
//     }
//   }
//   return x;
//
// }


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
# for example,
# try(inf_mat_match(c(0,2,4), c(2,5, 7), rnorm(10), 100))
*/
