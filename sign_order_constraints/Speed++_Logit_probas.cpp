// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
using namespace arma; // use the Armadillo library for matrix computations
using namespace Rcpp;

//Market Share Computations "Logit probas"
//[[Rcpp::export]]
mat probabilities_Logit_cpp(mat const& beta_draw, mat const& X){
  int draws = beta_draw.n_rows;
  int nplay = X.n_rows;
  mat Xbeta = zeros(draws,nplay);
  //Xbeta
  Xbeta = beta_draw * trans(X);
  mat market_share_draws = zeros(draws,nplay);
  //Compute stabilized choice probs now...
  for(int r = 0; r<draws; r++){
    vec Xbeta_temp = zeros(nplay); 
    Xbeta_temp = arma::conv_to< vec >::from(Xbeta.row(r)); //first row of matrix is NOT interpreted as vec
    vec Xbeta_temp_stab = zeros(nplay);
    vec max_temp = zeros(nplay);
    max_temp.fill(max(Xbeta_temp));
    Xbeta_temp_stab = Xbeta_temp - max_temp;
    vec numerator = exp(Xbeta_temp_stab);
    vec denominator = zeros(nplay);
    denominator.fill(sum(exp(Xbeta_temp_stab)));
    market_share_draws(r,span()) = trans(numerator/denominator);
  }
  return market_share_draws;
}

