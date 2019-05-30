// [[Rcpp::depends(RcppArmadillo)]]
#include "bayesm_LogN&normal.h"
#include "RcppArmadillo.h"
#include <stdio.h>
#include <time.h>
using namespace arma; // use the Armadillo library for matrix computations
using namespace Rcpp;


//[[Rcpp::export]]
vec startobeta(vec const& betastar){
  
  // by Max Pachali & Thomas Otter (2016-06-11)
  
  // converts normally distributed betastars to constrained betas
  
  int nvar = betastar.size();
  
  vec beta = zeros<vec>(nvar);
  
  // Price
  beta[0] = -exp(betastar[0]);
  // Cage (unconstrained)
  beta[14] = betastar[14];
  // Barn
  beta[1] = beta[14] + exp(betastar[1]);
  // Free Range
  beta[2] = beta[1] + exp(betastar[2]);
  // Organic
  beta[3] = beta[2] + exp(betastar[3]);
  // 6 Size
  beta[4] = -exp(betastar[4]);
  // Cage (unconstrained) Eastern
  beta[24] = betastar[24];
  // Cage (unconstrained) Xmas
  beta[25] = betastar[25];
  // Barn Eastern
  beta[5] = beta[24] + exp(betastar[5]);
  // Barn Xmas
  beta[6] = beta[25] + exp(betastar[6]);
  // Free-range Eastern
  beta[7] = beta[5] + exp(betastar[7]);
  // Free-range Xmas
  beta[8] = beta[6] + exp(betastar[8]);
  // Organic Eastern
  beta[9] = beta[7] + exp(betastar[9]);
  // Organic Xmas
  beta[10] = beta[8] + exp(betastar[10]);
  // Barn post ban
  beta[11] = beta[14] + exp(betastar[11]);
  // Free-range post ban
  beta[12] = beta[11] + exp(betastar[12]);
  // Organic post ban
  beta[13] = beta[12] + exp(betastar[13]);
  // Store parameters (unconstrained)
  beta.subvec(15,23) = betastar.subvec(15,23);
  return beta;
}

// [[Rcpp::export]]
mat subfrommat(mat const& X, uvec const& target){
  // Max Pachali 23/06/2016 
  // Function to select rows from matrix to create submatrix
  int num_cols = X.n_cols;
  int rows_target = target.size();
  mat submat = zeros(rows_target,num_cols);
  for(int i=0; i<rows_target; i++){
    submat(i,span::all) = X(target(i),span::all);
  }
  return(submat);
} 

// [[Rcpp::export]]
vec subfromvec(vec const& x, uvec const& target){
  // Max Pachali 23/06/2016 
  // Function to select elements from vector to create subvector
  int size_target = target.size();
  vec subvector = zeros(size_target);
  for(int i=0; i<size_target; i++){
    subvector(i) = x(target(i));
  }
  return(subvector);
} 


// [[Rcpp::export]]
double llmnl_varyCS(vec const& betastar, vec const& y, mat const& X, vec const& choice, vec const& p){
  
  // Wayne Taylor 9/7/2014
  
  // Modified by Max Pachali & Marco  Kotschedoff  (2016)
  
  // Evaluates log-likelihood for the multinomial logit model with varying choice sets
  
  // choice is a vector of zero & ones indicating which alternative is chosen in each choice set
  
  // p is a vector indicating the size of the choice set associated with each observation
  
  //int n = y.size(); // same as before
  
  vec choices = unique(p); // identify unique sizes of choice sets
  
  int n_choices = choices.size(); // number of unique sizes
  
  vec beta = startobeta(betastar);
  
  mat Xbeta = X*beta; // same as before: X%*%beta: should result in a (n*j x 1) vector
  vec llmnl_part = zeros(n_choices); 
  for(int k = 0; k<n_choices; k++){
    uvec same_size = find(p == choices(k)); // identify observation belonging to a choice set of same length
    int n = same_size.size()/choices(k); // number of choice sets being of same length of choices[k]
    int j = choices(k);
    mat Xbeta_new = subfrommat(Xbeta,same_size);
    vec xby = zeros<vec>(n); 
    vec denom = zeros<vec>(n);
    vec sub_choice = subfromvec(choice,same_size);
    for(int i = 0; i<n;i++){      
      for(int l = 0; l<j;l++){
        denom(i) = denom(i) + exp(Xbeta_new(i*j+l));
      }
      uvec y_mod = find(sub_choice.subvec(i*j,i*j+j-1)==1);
      xby(i) = Xbeta_new(i*j+y_mod(0));
    }
    llmnl_part(k) = sum(xby - log(denom));
  }
  return(sum(llmnl_part));
}


//Used in rhierMnlDP and rhierMnlRwMixture------------------------------------------------------------------------
mnlMetropOnceOut mnlMetropOnce_mod(vec const& y, mat const& X, vec const& oldbeta, 
                                   double oldll,vec s, mat const& incroot, 
                                   vec const& betabar, mat const& rootpi, vec const& choice, vec const& p){ 
  // Wayne Taylor 10/01/2014
  
  // function to execute rw metropolis for the MNL
  // y is n vector with element = 1,...,j indicating which alt chosen
  // X is nj x k matrix of xvalues for each of j alt on each of n occasions
  // RW increments are N(0,s^2*t(inc.root)%*%inc.root)
  // prior on beta is N(betabar,Sigma)  Sigma^-1=rootpi*t(rootpi)
  //  inc.root, rootpi are upper triangular
  //  this means that we are using the UL decomp of Sigma^-1 for prior 
  // oldbeta is the current
  
  // Likelihood modified by Max Pachali & Thomas Otter to include sign & order constraints
  
  mnlMetropOnceOut metropout_struct;
  
  double unif;
  vec betadraw, alphaminv;
  
  int stay = 0;
  int nvar = incroot.n_cols;
  mat s_incroot(nvar,nvar);
  for(int i=0; i<nvar; i++){
    s_incroot(span(),i) = s%incroot(span(),i);
  }
  vec betac = oldbeta + trans(s_incroot)*as<vec>(rnorm(X.n_cols));
  double cll = llmnl_varyCS(betac,y,X,choice,p);
  double clpost = cll+lndMvn(betac,betabar,rootpi);
  double ldiff = clpost-oldll-lndMvn(oldbeta,betabar,rootpi);
  alphaminv << 1 << exp(ldiff);
  double alpha = min(alphaminv);
  
  if(alpha < 1) {
    unif = runif(1)[0]; //runif returns a NumericVector, so using [0] allows for conversion to double
  } else { 
    unif=0;}
  if (unif <= alpha) {
    betadraw = betac;
    oldll = cll;
  } else {
    betadraw = oldbeta;
    stay = 1;
  }
  
  metropout_struct.betadraw = betadraw;
  metropout_struct.stay = stay;  
  metropout_struct.oldll = oldll;
  
  return (metropout_struct);
}


//[[Rcpp::export]]
List rhierMnlRwMixture_rcpp_loop_mod(List const& lgtdata, int nvar_c, mat const& Z, vec const& deltabar, mat const& Ad, mat const& mustarbarc, mat const& Amu,
                                     int const& nu, mat const& V, mat const& Gammabar, mat const& A_Gamma, int const& nu_Sigma, mat const& V_Sigma,
                                     vec s, int R, int keep, int nprint, bool drawdelta, mat olddelta,  vec const& a, vec oldprob, mat oldbetas, vec ind,
                                     mat const& betas_tuning, vec const& betastarpooled, vec ind_stay_temp, mat ind_stay, int ncomp){
  
  // Wayne Taylor 10/01/2014
  
  // Modified by Max Pachali & Thomas Otter to include sign & order constraints (log-normal prior)
  
  int nlgt = lgtdata.size();
  int nvar = oldbetas.n_cols;
  int nz = Z.n_cols;
  
  mat rootpi, betabar, ucholinv, incroot;
  int mkeep;
  mnlMetropOnceOut metropout_struct;
  List lgtdatai, nmix;
  
  // convert List to std::vector of struct
  std::vector<moments> lgtdata_vector;
  moments lgtdatai_struct;
  for (int lgt = 0; lgt<nlgt; lgt++){
    lgtdatai = lgtdata[lgt];
    
    lgtdatai_struct.y = as<vec>(lgtdatai["y"]);
    lgtdatai_struct.X = as<mat>(lgtdatai["X"]);
    lgtdatai_struct.hess = as<mat>(lgtdatai["hess"]);
    lgtdatai_struct.choice = as<vec>(lgtdatai["choice"]);
    lgtdatai_struct.p = as<vec>(lgtdatai["p"]);
    lgtdata_vector.push_back(lgtdatai_struct);    
  }
  
  // allocate space for draws
  vec oldll = zeros<vec>(nlgt);
  cube betadraw(nlgt, nvar, R/keep);
  mat probdraw(R/keep, oldprob.size());
  vec loglike(R/keep);
  mat Deltadraw(1,1); if(drawdelta) Deltadraw.zeros(R/keep, nz*nvar);//enlarge Deltadraw only if the space is required
  List compdraw(R/keep);
  // initialize rootic
  List prev_rootic(ncomp);
  mat rootic = eye(nvar_c,nvar_c);
  for(int f=0; f<ncomp; f++){
    prev_rootic(f) = List::create(
      Named("rootic") =  rootic);
  }
  
  if (nprint>0) startMcmcTimer();
  
  for (int rep = 0; rep<R; rep++){
    
    // first draw comps,ind,p | {beta_i}, delta
    // ind,p need initialization comps is drawn first in sub-Gibbs
    List mgout;
    if(drawdelta) {
      olddelta.reshape(nvar,nz);
      mgout = rmixGibbs_mod(oldbetas-Z*trans(olddelta),nvar_c,mustarbarc,Amu,nu,V,Gammabar,A_Gamma,nu_Sigma,V_Sigma,a,oldprob,ind,prev_rootic);
    } else {
      mgout = rmixGibbs_mod(oldbetas,nvar_c,mustarbarc,Amu,nu,V,Gammabar,A_Gamma,nu_Sigma,V_Sigma,a,oldprob,ind,prev_rootic);
    }
    
    List oldcomp = mgout["comps"];
    oldprob = as<vec>(mgout["p"]); //conversion from Rcpp to Armadillo requires explict declaration of variable type using as<>
    ind = as<vec>(mgout["z"]);
    
    //Update prev_rootic
    //List prev_rootic(ncomp);
    for(int f=0; f<ncomp; f++){
      List oldcomptemp = oldcomp[f];
      prev_rootic(f) = List::create(
        Named("rootic") =  as<mat>(oldcomptemp["rootic"]));
    }
    
    //now draw delta | {beta_i}, ind, comps
    if(drawdelta) olddelta = drawDelta(Z,oldbetas,ind,oldcomp,deltabar,Ad);
    
    //loop over all LGT equations drawing beta_i | ind[i],z[i,],mu[ind[i]],rooti[ind[i]]
    for(int lgt = 0; lgt<nlgt; lgt++){
      List oldcomplgt = oldcomp[ind[lgt]-1];
      rootpi = as<mat>(oldcomplgt[1]);
      
      //note: beta_i = Delta*z_i + u_i  Delta is nvar x nz
      if(drawdelta){
        olddelta.reshape(nvar,nz);
        betabar = as<vec>(oldcomplgt[0])+olddelta*vectorise(Z(lgt,span::all));
      } else {
        betabar = as<vec>(oldcomplgt[0]);
      }
      
      if (rep == 0) oldll[lgt] = llmnl_varyCS(vectorise(oldbetas(lgt,span::all)),lgtdata_vector[lgt].y,lgtdata_vector[lgt].X,lgtdata_vector[lgt].choice,lgtdata_vector[lgt].p);
      
      //compute inc.root
      //Note that lgtdata_vector[lgt].hess now represents the approximated Hessian (Delta method) 
      ucholinv = solve(trimatu(chol(lgtdata_vector[lgt].hess+rootpi*trans(rootpi))), eye(nvar,nvar)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
      incroot = chol(ucholinv*trans(ucholinv));
      
      metropout_struct = mnlMetropOnce_mod(lgtdata_vector[lgt].y,lgtdata_vector[lgt].X,vectorise(oldbetas(lgt,span::all)),
                                           oldll[lgt],s,incroot,betabar,rootpi,lgtdata_vector[lgt].choice,lgtdata_vector[lgt].p);
      
      oldbetas(lgt,span::all) = trans(metropout_struct.betadraw);
      oldll[lgt] = metropout_struct.oldll;
      ind_stay_temp[lgt] = metropout_struct.stay;
    }
    
    //print time to completion and draw # every nprint'th draw
    if (nprint>0) if ((rep+1)%nprint==0) infoMcmcTimer(rep, R);
    
    if((rep+1)%keep==0){
      mkeep = (rep+1)/keep;
      betadraw.slice(mkeep-1) = oldbetas;
      probdraw(mkeep-1, span::all) = trans(oldprob);
      loglike[mkeep-1] = sum(oldll);
      ind_stay(mkeep-1,span::all) = trans(ind_stay_temp);
      if(drawdelta) Deltadraw(mkeep-1, span::all) = trans(vectorise(olddelta));
      compdraw[mkeep-1] = oldcomp;
    }
  }
  
  if (nprint>0) endMcmcTimer();
  
  nmix = List::create(Named("probdraw") = probdraw,
                      Named("zdraw") = R_NilValue, //sets the value to NULL in R
                      Named("compdraw") = compdraw);
  
  if(drawdelta){
    return(List::create(
        Named("Deltadraw") = Deltadraw,
        Named("betadraw") = betadraw,
        Named("nmix") = nmix,
        Named("loglike") = loglike));  
  } else {
    return(List::create(
        Named("betadraw") = betadraw,
        Named("nmix") = nmix,
        Named("loglike") = loglike,
        Named("tuning") = List::create(Named("pooledstars") = betastarpooled, Named("individualstars")= betas_tuning),
              Named("rejection") = ind_stay));
  }
}


// [[Rcpp::export]]
List rmultireg(mat const& Y, mat const& X, mat const& Bbar, mat const& A, int nu, mat const& V) {
  
  // Keunwoo Kim 09/09/2014
  
  // Purpose: draw from posterior for Multivariate Regression Model with natural conjugate prior
  
  // Arguments:
  //  Y is n x m matrix
  //  X is n x k
  //  Bbar is the prior mean of regression coefficients  (k x m)
  //  A is prior precision matrix
  //  nu, V are parameters for prior on Sigma
  
  // Output: list of B, Sigma draws of matrix of coefficients and Sigma matrix
  
  // Model: 
  //  Y=XB+U  cov(u_i) = Sigma
  //  B is k x m matrix of coefficients
  
  // Prior:  
  //  beta|Sigma  ~ N(betabar,Sigma (x) A^-1)
  //  betabar=vec(Bbar)
  //  beta = vec(B) 
  //  Sigma ~ IW(nu,V) or Sigma^-1 ~ W(nu, V^-1)
  
  int n = Y.n_rows;
  int m = Y.n_cols;
  int k = X.n_cols;
  
  //first draw Sigma
  mat RA = chol(A);
  mat W = join_cols(X, RA); //analogous to rbind() in R
  mat Z = join_cols(Y, RA*Bbar);
  // note:  Y,X,A,Bbar must be matrices!
  mat IR = solve(trimatu(chol(trans(W)*W)), eye(k,k)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
  // W'W = R'R  &  (W'W)^-1 = IRIR'  -- this is the UL decomp!
  mat Btilde = (IR*trans(IR)) * (trans(W)*Z);
  // IRIR'(W'Z) = (X'X+A)^-1(X'Y + ABbar)
  mat E = Z-W*Btilde;
  mat S = trans(E)*E;
  // E'E
  
  // compute the inverse of V+S
  mat ucholinv = solve(trimatu(chol(V+S)), eye(m,m));
  mat VSinv = ucholinv*trans(ucholinv);
  
  List rwout = rwishart(nu+n, VSinv);
  
  // now draw B given Sigma
  //   note beta ~ N(vec(Btilde),Sigma (x) Covxxa)
  //       Cov=(X'X + A)^-1  = IR t(IR)  
  //       Sigma=CICI'    
  //       therefore, cov(beta)= Omega = CICI' (x) IR IR' = (CI (x) IR) (CI (x) IR)'
  //  so to draw beta we do beta= vec(Btilde) +(CI (x) IR)vec(Z_mk)  
  //  		Z_mk is m x k matrix of N(0,1)
  //	since vec(ABC) = (C' (x) A)vec(B), we have 
  //		B = Btilde + IR Z_mk CI'
  
  mat CI = rwout["CI"]; //there is no need to use as<mat>(rwout["CI"]) since CI is being initiated as a mat in the same line
  mat draw = mat(rnorm(k*m));
  draw.reshape(k,m);
  mat B = Btilde + IR*draw*trans(CI);
  
  return List::create(
    Named("B") = B, 
    Named("Sigma") = rwout["IW"]);
}


// [[Rcpp::export]]
List rwishart(int const& nu, mat const& V){
  
  // Wayne Taylor 4/7/2015
  
  // Function to draw from Wishart (nu,V) and IW
  
  // W ~ W(nu,V)
  // E[W]=nuV
  
  // WI=W^-1
  // E[WI]=V^-1/(nu-m-1)
  
  // T has sqrt chisqs on diagonal and normals below diagonal
  int m = V.n_rows;
  mat T = zeros(m,m);
  
  for(int i = 0; i < m; i++) {
    T(i,i) = sqrt(rchisq(1,nu-i)[0]); //rchisq returns a vectorized object, so using [0] allows for the conversion to double
  }
  
  for(int j = 0; j < m; j++) {  
    for(int i = j+1; i < m; i++) {    
      T(i,j) = rnorm(1)[0]; //rnorm returns a NumericVector, so using [0] allows for conversion to double
    }}
  
  mat C = trans(T)*chol(V);
  mat CI = solve(trimatu(C),eye(m,m)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
  
  // C is the upper triangular root of Wishart therefore, W=C'C
  // this is the LU decomposition Inv(W) = CICI' Note: this is
  // the UL decomp not LU!
  
  // W is Wishart draw, IW is W^-1
  
  return List::create(
    Named("W") = trans(C) * C,
    Named("IW") = CI * trans(CI),
    Named("C") = C,
    Named("CI") = CI);
}


// rmixGibbs
List drawCompsFromLabels(mat const& y, int nvar_c, mat const& mustarbarc, mat const& Amu, int nu, mat const& V, mat const& Gammabar, mat const& A_Gamma,  
                          int nu_Sigma, mat const& V_Sigma, int ncomp, vec const& z, List prev_rootic){
  
  // Wayne Taylor 3/18/2015
  // Function to draw the components based on the z labels
  
  // Modified by Max Pachali 
  // Separate update of hierarchical prior on beta^c and beta^uc
  
  // Gammabar is (nvar_c x nvar - nvar_c) & A_Gamma is (nvar_c x nvar_c) 
  // & V_Sigma is (nvar - nvar_c x nvar - nvar_c)
  // mustarbarc is (nvar_c x 1) & Vbetastarc is (nvar_c x nvar_c)
  
  List comps(ncomp);
  
  int n = z.n_rows;
  int nvar = y.n_cols;
  int nvar_uc = nvar - nvar_c;
  vec nobincomp = zeros<vec>(ncomp);
  
  //Determine the number of observations in each component
  for(int i = 0; i<n; i++) {
    nobincomp[z[i]-1]++; //Note z starts at 1, not 0
  }
  
  //Draw comps
  //Note that y is a matrix of dimension (Nxnvar)
  for(int k = 0; k<ncomp; k++){
    
      // If there are observations in this component, draw from the posterior
      if(nobincomp[k] > 0) {
      mat yk = y.rows(find(z==(k+1))); //Note k starts at 0 and z starts at 1
      // Divide the y's into constrained & unconstrained observations (requires ordered design matrix into constrained & unconstrained coefficients)
      mat yk_c = yk(span::all,span(0,nvar_c-1));
      mat consta = ones(nobincomp[k],1);
      mat yk_c_const = join_rows(consta,yk_c);
      mat yk_uc = yk(span::all,span(nvar_c,nvar-1));
      
      ////////////////////////////////////////////////////////////////////////////////
      // Update Gamma & Sigma by regressing yk_c on yk_uc (natural conjugate draw) //
      //////////////////////////////////////////////////////////////////////////////
      List temp1 = rmultireg(yk_uc, yk_c_const, Gammabar, A_Gamma, nu_Sigma, V_Sigma);
      mat Sigma = as<mat>(temp1["Sigma"]); //conversion from Rcpp to Armadillo requires explict declaration of variable type using as<>
      mat Gamma = as<mat>(temp1["B"]); //Note that the first row of Gamma are intercept terms
      
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Update mustarc & Vbetastarc by regressing a vector of ones on betastardraw_c with CONDITIONALLY conjugate priors //
      // Use part from DrawDelta-function to update mustarc conditional on Vbetastarc /////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
      //mustarc given Vbetastarc first
      mat Xk = ones(nobincomp[k], 1);
      //Xk.reshape(nobincomp[k],1);
      mat xtx = zeros<mat>(nvar_c,nvar_c);
      mat xty = zeros<mat>(nvar_c,1);
      //Get rootic
      List rootictemp = prev_rootic(k);
      mat rootic = as<mat>(rootictemp["rootic"]);
      mat Vbetastarc_inv = rootic*trans(rootic);
      xtx = kron((trans(Xk)*Xk),Vbetastarc_inv);
      xty = Vbetastarc_inv * (trans(yk_c)*Xk);
      // compute the inverse of xtx+A
      mat ucholinv = solve(trimatu(chol(xtx+Amu)), eye(nvar_c,nvar_c)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
      mat Ainv_draw = ucholinv*trans(ucholinv);
      // Update mustarc now
      mat mustarc_draw = Ainv_draw*(xty + Amu*mustarbarc) + trans(chol(Ainv_draw))*as<vec>(rnorm(mustarbarc.size()));
      
      // Vbetastarc given fresh mustarc draw
      int nu_update = nu + nobincomp[k];
      mat mustarc_drawN = Xk*trans(mustarc_draw);
      mat E = yk_c-mustarc_drawN;
      mat ykmu = trans(E)*E;
      mat iVupdate_r = solve(trimatu(chol(V+ykmu)), eye(V.n_rows,V.n_cols));
      mat iVupdate = iVupdate_r*trans(iVupdate_r);
      // Draw Vbetastarc now
      List rwout = rwishart(nu_update, iVupdate);
      mat Vbetastarc_draw = as<mat>(rwout["IW"]);
      //Update rootic
      mat rootic_up = solve(trimatu(chol(Vbetastarc_draw)),eye(Vbetastarc_draw.n_rows,Vbetastarc_draw.n_cols));
      
      //Create aggregate Vbetastar from ingredients
      mat Vbetastar = zeros(nvar,nvar);
      mat Vbetastarcgamma = Vbetastarc_draw*Gamma(span(1,nvar_c),span::all);
      mat Vbetastaruc = (trans(Gamma(span(1,nvar_c),span::all))*Vbetastarc_draw)*Gamma(span(1,nvar_c),span::all) + Sigma;
      Vbetastar(span(0,nvar_c-1),span(0,nvar_c-1)) = Vbetastarc_draw;
      Vbetastar(span(0,nvar_c-1),span(nvar_c,nvar-1)) = Vbetastarcgamma;
      Vbetastar(span(nvar_c,nvar-1),span(0,nvar_c-1)) = trans(Vbetastarcgamma);
      Vbetastar(span(nvar_c,nvar-1),span(nvar_c,nvar-1)) = Vbetastaruc;
      mat rooti = solve(trimatu(chol(Vbetastar)),eye(Vbetastar.n_rows,Vbetastar.n_cols));
      //Create aggregate mu from ingredients
      mat mustar = zeros(nvar,1);
      mat mustarcconst_draw = ones(nvar_c+1,1);
      mustarcconst_draw(span(1,nvar_c),0) = mustarc_draw;
      mat Gammamustar = trans(Gamma)*mustarcconst_draw;
      Gammamustar.reshape(nvar_uc,1);
      mustar(span(0,nvar_c-1),0) = mustarc_draw;
      mustar(span(nvar_c,nvar-1),0) = Gammamustar;

      comps(k) = List::create(
        Named("mustar") = NumericVector(mustar.begin(),mustar.end()),
        Named("rooti") = rooti,
        Named("rootic") = rootic_up
      );
      
    } else {
      // If there are no obervations in this component, draw from the prior
      
      // Update Sigma & Gamma 
      mat S = solve(trimatu(chol(V_Sigma)),eye(V_Sigma.n_rows,V_Sigma.n_cols));
      S = S * trans(S); 
      List rw = rwishart(nu_Sigma, S);
      mat Sigma = as<mat>(rw["IW"]);
      mat CI = as<mat>(rw["CI"]);
      // Gamma
      int dim = Gammabar.n_rows;
      int m = Gammabar.n_cols;
      mat RA = chol(A_Gamma);
      mat W = RA;
      mat IR = solve(trimatu(chol(trans(W)*W)), eye(dim,dim)); 
      mat draw = mat(rnorm(dim*m));
      draw.reshape(dim,m);
      mat Gamma = Gammabar + IR*draw*trans(CI);
      
      // Update mustarc & Vbetastarc
      mat ucholinv = solve(trimatu(chol(Amu)), eye(Amu.n_rows,Amu.n_cols)); 
      mat Ainv_draw = ucholinv*trans(ucholinv);
      mat mustarc_draw = mustarbarc + trans(chol(Ainv_draw))*as<vec>(rnorm(mustarbarc.size()));
      mustarc_draw.reshape(nvar_c,1);
      
      // Vbetastarc
      mat V_inv = solve(trimatu(chol(V)), eye(V.n_rows,V.n_cols));
      List rwout = rwishart(nu, V_inv);
      mat Vbetastarc_draw = as<mat>(rwout["IW"]);
      //Update rootic
      mat rootic_up = solve(trimatu(chol(Vbetastarc_draw)),eye(Vbetastarc_draw.n_rows,Vbetastarc_draw.n_cols));
      
      //Create aggregate Vbetastar from ingredients
      mat Vbetastar = zeros(nvar,nvar);
      mat Vbetastarcgamma = Vbetastarc_draw*Gamma(span(1,nvar_c),span::all);
      mat Vbetastaruc = (trans(Gamma(span(1,nvar_c),span::all))*Vbetastarc_draw)*Gamma(span(1,nvar_c),span::all) + Sigma;
      Vbetastar(span(0,nvar_c-1),span(0,nvar_c-1)) = Vbetastarc_draw;
      Vbetastar(span(0,nvar_c-1),span(nvar_c,nvar-1)) = Vbetastarcgamma;
      Vbetastar(span(nvar_c,nvar-1),span(0,nvar_c-1)) = trans(Vbetastarcgamma);
      Vbetastar(span(nvar_c,nvar-1),span(nvar_c,nvar-1)) = Vbetastaruc;
      mat rooti = solve(trimatu(chol(Vbetastar)),eye(Vbetastar.n_rows,Vbetastar.n_cols));
 
      //Create aggregate mu from ingredients
      mat mustar = zeros(nvar,1);
      mat mustarcconst_draw = ones(nvar_c+1,1);
      mustarcconst_draw(span(0,nvar_c-1),0) = mustarc_draw;
      mat Gammamustar = trans(Gamma)*mustarcconst_draw;
      Gammamustar.reshape(nvar_uc,1);
      mustar(span(0,nvar_c-1),0) = mustarc_draw;
      mustar(span(nvar_c,nvar-1),0) = Gammamustar;
      
      comps(k) = List::create(
      Named("mustar") = NumericVector(mustar.begin(),mustar.end()),
      Named("rooti") = rooti,
      Named("rootic") = rootic_up
      );
    } 
  }
  
  return(comps);
}

vec drawLabelsFromComps(mat const& y, vec const& p, List comps) {
  
  // Wayne Taylor 3/18/2015
  
  // Function to determine which label is associated with each y value
  
  double logprod;
  vec mu, u;
  mat rooti;
  List compsk;
  
  int n = y.n_rows;
  vec res = zeros<vec>(n);
  int ncomp  = comps.size();
  mat prob(n,ncomp);
  
  for(int k = 0; k<ncomp; k++) {
    compsk = comps[k];
    mu = as<vec>(compsk["mustar"]); //conversion from Rcpp to Armadillo requires explict declaration of variable type using as<>
    rooti = as<mat>(compsk["rooti"]);
    
    //Find log of MVN density using matrices
    logprod = log(prod(diagvec(rooti)));
    mat z(y);
    z.each_row() -= trans(mu); //subtracts mu from each row in z
    z = trans(rooti) * trans(z);
    z = -(y.n_cols/2.0) * log(2*M_PI) + logprod - .5 * sum(z % z, 0); // operator % performs element-wise multiplication
    
    prob.col(k) =  trans(z);
  }
  
  prob = exp(prob);
  prob.each_row() %= trans(p); //element-wise multiplication
  
  // Cumulatively add each row and take a uniform draw between 0 and the cumulative sum
  prob = cumsum(prob, 1);
  u = as<vec>(runif(n)) % prob.col(ncomp-1);
  
  // Evaluative each column of "prob" until the uniform draw is less than the cumulative value
  for(int i = 0; i<n; i++) {
    while(u[i] > prob(i, res[i]++));
  }
  
  return(res);
}

// [[Rcpp::export]]
vec rdirichlet(vec const& alpha){
  
  // Wayne Taylor 4/7/2015
  
  // Purpose:
  // draw from Dirichlet(alpha)
  
  int dim = alpha.size();
  vec y = zeros<vec>(dim);
  
  for(int i = 0; i<dim; i++) {    
    y[i] = rgamma(1,alpha[i])[0]; //rgamma returns a NumericVector, so adding [0] extracts the first element and treats it as type "double"
  }
  
  return(y/sum(y));
}


vec drawPFromLabels(vec const& a, vec const& z) {
  
  // Wayne Taylor 9/10/2014
  
  // Function to draw the probabilities based on the label proportions
  
  vec a2 = a;
  int n = z.n_rows;
  
  //Count number of observations in each component
  for(int i = 0; i<n; i++) a2[z[i]-1]++; //z starts at 1, not 0
  
  return rdirichlet(a2);
}

//[[Rcpp::export]]
List rmixGibbs_mod( mat const& y, int nvar_c, mat const& mustarbarc, mat const& Amu, int nu, mat const& V, mat const& Gammabar, mat const& A_Gamma,
                int nu_Sigma, mat const& V_Sigma, vec const& a, vec const& p,  vec const& z, List prev_rootic) {
  
  // Wayne Taylor 9/10/2014
  
  /*
  // Revision History: R. McCulloch 11/04 P. Rossi 3/05 put in
  // backsolve and improved documentation
  // 
  // purpose: do gibbs sampling inference for a mixture of
  // multivariate normals
  // 
  // arguments: y: data, rows are observations, assumed to be iid
  // draws from normal mixture Bbar,A,nu,V: common prior for mean
  // and variance of each normal component
  // 
  // note: Bbar should be a matrix. usually with only one row
  // 
  // beta ~ N(betabar,Sigma (x) A^-1) betabar=vec(Bbar) Sigma ~
  // IW(nu,V) or Sigma^-1 ~ W(nu, V^-1) note: if you want Sigma ~
  // A, use nu big and rwishart(nu,nu(A)^{-1})$IW a: Dirichlet
  // parameters for prior on p p: prior probabilities of normal
  // components z: components indentities for each observation
  // (vector of intergers each in {1,2,...number of components})
  // comps: list, each member is a list comp with ith normal
  // component ~N(comp[[1]],Sigma), Sigma = t(R)%*%R, R^{-1} =
  // comp[[2]] Output: list with elements [[1]=$p, [[2]]=$z, and
  // [[3]]=$comps, with the updated values
  
  */
  
  List comps = drawCompsFromLabels(y, nvar_c, mustarbarc, Amu, nu, V, Gammabar, A_Gamma, nu_Sigma, V_Sigma, a.size(), z, prev_rootic);

  vec z2 = drawLabelsFromComps(y, p, comps);
  
  vec p2 = drawPFromLabels(a, z2);
  
  return List::create(
    Named("p") = p2,
    Named("z") = z2,
    Named("comps") = comps);
}

// [[Rcpp::export]]
double lndMvn(vec const& x, vec const& mu, mat const& rooti){
  
  //Wayne Taylor 9/7/2014
  
  // function to evaluate log of MV Normal density with  mean mu, var Sigma
  // Sigma=t(root)%*%root   (root is upper tri cholesky root)
  // Sigma^-1=rooti%*%t(rooti)   
  // rooti is in the inverse of upper triangular chol root of sigma
  //          note: this is the UL decomp of sigmai not LU!
  //                Sigma=root'root   root=inv(rooti)
  
  vec z = vectorise(trans(rooti)*(x-mu));
  
  return((-(x.size()/2.0)*log(2*M_PI) -.5*(trans(z)*z) + sum(log(diagvec(rooti))))[0]);
}



//Used in rhierLinearModel, rhierLinearMixture and rhierMnlRWMixture------------------------------------------------------
// [[Rcpp::export]]
mat drawDelta(mat const& x,mat const& y,vec const& z,List const& comps,vec const& deltabar,mat const& Ad){
  
  // Wayne Taylor 10/01/2014
  
  // delta = vec(D)
  //  given z and comps (z[i] gives component indicator for the ith observation, 
  //   comps is a list of mu and rooti)
  // y is n x p
  // x is n x k
  // y = xD' + U , rows of U are indep with covs Sigma_i given by z and comps
  
  int p = y.n_cols;
  int k = x.n_cols;
  int ncomp  = comps.length();
  mat xtx = zeros<mat>(k*p,k*p);
  mat xty = zeros<mat>(p,k); //this is the unvecced version, reshaped after the sum
  
  //Create the index vectors, the colAll vectors are equal to span::all but with uvecs (as required by .submat)
  uvec colAlly(p), colAllx(k);
  for(int i = 0; i<p; i++) colAlly(i) = i;
  for(int i = 0; i<k; i++) colAllx(i) = i;
  
  //Loop through the components
  for(int compi = 0; compi<ncomp; compi++){
    
    //Create an index vector ind, to be used like y[ind,]
    uvec ind = find(z == (compi+1));
    
    //If there are observations in this component
    if(ind.size()>0){
      mat yi = y.submat(ind,colAlly);
      mat xi = x.submat(ind,colAllx);
      
      List compsi = comps[compi];
      rowvec mui = as<rowvec>(compsi[0]); //conversion from Rcpp to Armadillo requires explict declaration of variable type using as<>
      mat rootii = trimatu(as<mat>(compsi[1])); //trimatu interprets the matrix as upper triangular
      yi.each_row() -= mui; //subtracts mui from each row of yi
      mat sigi = rootii*trans(rootii);
      xtx = xtx + kron(trans(xi)*xi,sigi);
      xty = xty + (sigi * (trans(yi)*xi));
    }
  }
  xty.reshape(xty.n_rows*xty.n_cols,1);
  
  //vec(t(D)) ~ N(V^{-1}(xty + Ad*deltabar),V^{-1}) where V = (xtx+Ad)
  // compute the inverse of xtx+Ad
  mat ucholinv = solve(trimatu(chol(xtx+Ad)), eye(k*p,k*p)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
  mat Vinv = ucholinv*trans(ucholinv);
  
  return(Vinv*(xty+Ad*deltabar) + trans(chol(Vinv))*as<vec>(rnorm(deltabar.size())));
}


//The functions below are used to print the output from MCMC draws for many of the bayesm functions
time_t itime;
char buf[100];

//[[Rcpp::export]]
void startMcmcTimer() {
  itime = time(NULL);
  Rcout << " MCMC Iteration (est time to end - min) \n";
}

//[[Rcpp::export]]
void infoMcmcTimer(int rep, int R) {
  time_t ctime = time(NULL);    
  char buf[32];
  
  double timetoend = difftime(ctime, itime) / 60.0 * (R - rep - 1) / (rep+1);
  sprintf(buf, " %d (%.1f)\n", rep+1, timetoend);
  Rcout <<  buf;
}

//[[Rcpp::export]]
void endMcmcTimer() {
  time_t ctime = time(NULL);
  char buf[32];
  
  sprintf(buf, " Total Time Elapsed: %.2f \n", difftime(ctime, itime) / 60.0);     
  Rcout << buf;
  
  itime = 0;
}

