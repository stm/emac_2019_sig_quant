rhierMnlRwMixture_SR <- function (Data, Prior, Mcmc, nvar_c) 
{
  fsh=function() 
  {
    # 
    # P. Rossi
    # revision history: 3/27/05
    #
    # Purpose:
    #  function to flush console (needed only under windows)
    #
    if (Sys.info()[1] == "Windows") flush.console()
    return()
  }
  
  pandterm=function(message) {stop(message,call.=FALSE)}
  
  if (missing(Data)) {
    pandterm("Requires Data argument -- list of p,lgtdata, and (possibly) Z")
  }
  #if (is.null(Data$p)) {
  #pandterm("Requires Data element p (# chce alternatives)")
  #}
  #p = Data$p
  #if (is.null(Data$lgtdata)) {
  #pandterm("Requires Data element lgtdata (list of data for each unit)")
  #}
  lgtdata = Data$lgtdata
  nlgt = length(lgtdata)
  drawdelta = TRUE
  if (is.null(Data$Z)) {
    cat("Z not specified", fill = TRUE)
    #fsh()
    drawdelta = FALSE
  }
  else {
    if (nrow(Data$Z) != nlgt) {
      pandterm(paste("Nrow(Z) ", nrow(Z), "ne number logits ", 
                     nlgt))
    }
    else {
      Z = Data$Z
    }
  }
  if (drawdelta) {
    nz = ncol(Z)
    colmeans = apply(Z, 2, mean)
    if (sum(colmeans) > 1e-05) {
      pandterm(paste("Z does not appear to be de-meaned: colmeans= ", 
                     colmeans))
    }
  }
  ypooled = NULL
  Xpooled = NULL
  ppooled=NULL
  choicepooled=NULL
  
  if (!is.null(lgtdata[[1]]$X)) {
    oldncol = ncol(lgtdata[[1]]$X)
  }
  for (i in 1:nlgt) {
    if(is.null(lgtdata[[i]]$y)) {pandterm(paste("Requires element y of lgtdata[[",i,"]]"))}
    if(is.null(lgtdata[[i]]$X)) {pandterm(paste("Requires element X of lgtdata[[",i,"]]"))}
    ###NEW!### 
    if(is.null(lgtdata[[i]]$p)) {pandterm(paste("Requires element p of lgtdata[[",i,"]]"))}
    #  
    ypooled=c(ypooled,lgtdata[[i]]$y)
    ppooled=c(ppooled,lgtdata[[i]]$p)
    nrowX=nrow(lgtdata[[i]]$X)
    #if((nrowX/p) !=length(lgtdata[[i]]$y)) {pandterm(paste("nrow(X) ne p*length(yi); exception at unit",i))}
    ###NEW!###
    if(length(lgtdata[[i]]$p)!=length(lgtdata[[i]]$choice)) {pandterm(paste("lenght(pi) ne *length(yi); exception at unit",i))}
    newncol=ncol(lgtdata[[i]]$X)
    if(newncol != oldncol) {pandterm(paste("All X elements must have same # of cols; exception at unit",i))}
    Xpooled=rbind(Xpooled,lgtdata[[i]]$X)
    choicepooled=c(choicepooled,lgtdata[[i]]$choice)
    oldncol=newncol
  }
  nvar = ncol(Xpooled)
  nvar_uc = nvar-nvar_c
  levely = as.numeric(levels(as.factor(ypooled)))
  #if (length(levely) != p) {
  #pandterm(paste("y takes on ", length(levely), " values -- must be = p"))
  #}
  #bady = FALSE
  #for (i in 1:p) {
  #if (levely[i] != i) 
  #bady = TRUE
  #}
  cat("Table of Y values pooled over all units", fill = TRUE)
  print(table(ypooled))
  #if (bady) {
  #pandterm("Invalid Y")
  #}
  
  if (missing(Prior)) {
    pandterm("Requires Prior list argument (at least ncomp)")
  }
  if (is.null(Prior$ncomp)) {
    pandterm("Requires Prior element ncomp (num of mixture components)")
  }else{
    ncomp = Prior$ncomp
  }
  if (is.null(Prior$mustarbarc)) {
    mustarbarc = matrix(rep(0, nvar_c), nrow = nvar_c)
  } else {
    mustarbarc = matrix(Prior$mustarbarc, nrow = nvar_c)
  }
  if (is.null(Prior$Gammabar)) {
    Gammabar = matrix(0, nrow = (nvar_c+1), ncol = nvar_uc)
  }  else {
    Gammabar = matrix(Prior$Gammabar, nrow = (nvar_c+1), ncol = nvar_uc)
  }
  if (is.null(Prior$Amu)) {
    Amu = diag(0.1, nrow = nvar_c, ncol = nvar_c)
    ###change the default prior for Amu
  } else {
    Amu = matrix(Prior$Amu, nrow = nvar_c, ncol = nvar_c)
  }
  if (is.null(Prior$A_Gamma)) {
    A_Gamma = diag(0.01, nrow = (nvar_c+1), ncol = (nvar_c+1))
  } else {
    A_Gamma = matrix(Prior$A_Gamma, nrow = (nvar_c+1), ncol = (nvar_c+1))
  }
  if (is.null(Prior$nu)) {
    nu = nvar_c + 15
    ###Informative prior on restricted coefficients
  } else {
    nu = Prior$nu
  }
  if (is.null(Prior$nu_Sigma)) {
    nu_Sigma = nvar_uc + 5
  }  else {
    nu_Sigma = Prior$nu_Sigma
  }
  if (is.null(Prior$V)) {
    V = nu * diag(nvar_c)*0.5
    ###Allenby et al. (2014) 
  } else {
    V = Prior$V
  }
  if (is.null(Prior$V_Sigma)) {
    V_Sigma = nu_Sigma * diag(nvar_uc)
  } else {
    V_Sigma = Prior$V_Sigma
  }
  if (is.null(Prior$Ad) & drawdelta) {
    Ad = 0.01 * diag(nvar * nz)
  } else {
    Ad = Prior$Ad
  }
  if (drawdelta) {
    if (ncol(Ad) != nvar * nz | nrow(Ad) != nvar * nz) {
      pandterm("Ad must be nvar*nz x nvar*nz")
    }
  }
  if (is.null(Prior$deltabar) & drawdelta) {
    deltabar = rep(0, nz * nvar)
  } else {
    deltabar = Prior$deltabar
  }
  if (drawdelta) {
    if (length(deltabar) != nz * nvar) {
      pandterm("deltabar must be of length nvar*nz")
    }
  }
  if (is.null(Prior$a)) {
    a = rep(5, ncomp)
  } else {
    a = Prior$a
  }
  if (length(a) != ncomp) {
    pandterm("Requires dim(a)= ncomp (no of components)")
  }
  bada = FALSE
  for (i in 1:ncomp) {
    if (a[i] < 0) 
      bada = TRUE
  }
  if (bada) 
    pandterm("invalid values in a vector")
  if (missing(Mcmc)) {
    pandterm("Requires Mcmc list argument")
  } else {
    if (is.null(Mcmc$s)) {
      s = 2.38/sqrt(nvar)
    } else {
      s = Mcmc$s
    }
    if (is.null(Mcmc$w)) {
      w = 0.1
    } else {
      w = Mcmc$w
    }
    if (is.null(Mcmc$keep)) {
      keep = 1
    } else {
      keep = Mcmc$keep
    }
    if (is.null(Mcmc$R)) {
      pandterm("Requires R argument in Mcmc list")
    } else {
      R = Mcmc$R
    }
    if (is.null(Mcmc$nprint)) {
      nprint = 100
    } else {
      nprint = Mcmc$nprint
    }
    if (nprint < 0) {
      pandterm("nprint must be an integer greater than or equal to 0")
    }
  }
  cat(" ", fill = TRUE)
  cat("Starting MCMC Inference for Hierarchical Logit:", fill = TRUE)
  cat("   Normal Mixture with", ncomp, "components for first stage prior", 
      fill = TRUE)
  #cat(paste("  ", p, " alternatives; ", nvar, " variables in X"), 
      #fill = TRUE)
  cat(paste("   for ", nlgt, " cross-sectional units"), fill = TRUE)
  cat(" ", fill = TRUE)
  cat("Prior Parms: ", fill = TRUE)
  cat("nu =", nu, "nu_Sigma =", nu_Sigma,  fill = TRUE)
  cat("V ", fill = TRUE)
  print(V)
  cat("V_Sigma ", fill = TRUE)
  print(V_Sigma)
  cat("mustarbarc ", fill = TRUE)
  print(mustarbarc)
  cat("Amu ", fill = TRUE)
  print(Amu)
  cat("Gammabar ", fill = TRUE)
  print(Gammabar)
  cat("A_Gamma ", fill = TRUE)
  print(A_Gamma)
  cat("a ", fill = TRUE)
  print(a)
  if (drawdelta) {
    cat("deltabar", fill = TRUE)
    print(deltabar)
    cat("Ad", fill = TRUE)
    print(Ad)
  }
  cat(" ", fill = TRUE)
  cat("MCMC Parms: ", fill = TRUE)
  cat(paste("s=", round(s, 3), " w= ", w, " R= ", R, " keep= ", 
            keep, " nprint= ", nprint), fill = TRUE)
  cat("", fill = TRUE)
  oldbetas = matrix(double(nlgt * nvar), ncol = nvar)
  ind_stay = array(0,dim=c((R/keep),nlgt))
  ind_stay_temp = array(0,dim=c(nlgt,1))
  
  ###NEW!###: llmnl_varyCS instead of llmnl, also additional parameter choice and ppooled
  llmnlFract = function(beta,y,X,betapooled,rootH,w,wgt,choice,p){
    z=as.vector(rootH%*%(beta-betapooled))
    return((1-w)*llmnl_varyCS(beta,y,X,choice,p)+w*wgt*(-.5*(z%*%z)))
  }
  
  mnlHess_varyCS = function(betastar,y,X,p) {
    #   p.rossi 2004
    #   changed argument order 9/05
    #
    # Purpose: compute mnl -Expected[Hessian]  
    #
    # Arguments:
    #   beta is k vector of coefs
    #   y is n vector with element = 1,...,j indicating which alt chosen
    #   X is nj x k matrix of xvalues for each of j alt on each of n occasions
    #
    # Output:  -Hess evaluated at beta
    #
    #choices=levels(as.factor(p))
    #n_choices=length(choices) # that different legth of choice sets
    n = length(y)
    k=ncol(X)
    Hess=matrix(double(k*k),ncol=k)
    beta = startobeta(betastar)
    #loop over different choice set lengths
    #for (i in 1:n_choices) {
    #n=length(y[p[choice==1]==choices[i]])
    #X_new=X[p==choices[i],]
    #j=nrow(X_new)/n    
    #Xbeta=X_new%*%beta
    #Xbeta=matrix(Xbeta,byrow=T,ncol=j)
    #Xbeta=exp(Xbeta)
    #iota=c(rep(1,j))
    #denom=Xbeta%*%iota
    #Prob=Xbeta/as.vector(denom)
    cum_ind = 0
    for (m in 1:n) {
      if(m == 1){
        j = p[m]
      }else{
        j = p[(cum_ind+1)]
      }
      Xbeta = X[(cum_ind+1):(cum_ind+j),]%*%beta
      Xbeta = matrix(Xbeta,byrow=T,ncol=j)
      Xbeta = exp(Xbeta)
      iota = c(rep(1,j))
      denom = Xbeta%*%iota
      Prob = Xbeta/as.vector(denom)
      A = diag(as.vector(Prob))-outer(as.vector(Prob),as.vector(Prob))
      Xt = X[(cum_ind+1):(cum_ind+j),]
      Hess=Hess+crossprod(Xt,A)%*%Xt # updating the hessian
      cum_ind = cum_ind + j
    }
    #Hess = as.matrix(forceSymmetric(Hess))
    return(Hess)
  }
  
  # ###Modified function for hessian
  # mnlHess_mod = function (betastar, y, X) {
  #   beta = startobeta(betastar)  # call to compiled transformation routine
  #   n = length(y)
  #   j = nrow(X)/n
  #   k = ncol(X)
  #   Xbeta = X %*% beta
  #   Xbeta = matrix(Xbeta, byrow = T, ncol = j)
  #   Xbeta = exp(Xbeta)
  #   iota = c(rep(1, j))
  #   denom = Xbeta %*% iota
  #   Prob = Xbeta/as.vector(denom)
  #   Hess = matrix(double(k * k), ncol = k)
  #   for (i in 1:n) {
  #     p = as.vector(Prob[i, ])
  #     A = diag(p) - outer(p, p)
  #     Xt = X[(j * (i - 1) + 1):(j * i), ]
  #     Hess = Hess + crossprod(Xt, A) %*% Xt
  #   }
  #   return(Hess)
  # }
  # 
  # ###Compute gradient
  # mnlgrad_mod = function (betastar, y, X) {
  #   beta = startobeta(betastar)  # call to compiled transformation routine
  #   n = length(y)
  #   j = nrow(X)/n
  #   k = ncol(X)
  #   Xbeta = X %*% beta
  #   Xbeta = matrix(Xbeta, byrow = T, ncol = j)
  #   Xbeta = exp(Xbeta)
  #   iota = c(rep(1, j))
  #   denom = Xbeta %*% iota
  #   Prob = Xbeta/as.vector(denom)
  #   grad = t(array(0,dim=c(k,1)))
  #   for (i in 1:n) {
  #     p = as.vector(Prob[i, ]) ###for every choice occasion
  #     y_ind = rep(0,j); y_ind[y[i]] = 1 ###rewrite data
  #     Xt = X[(j * (i - 1) + 1):(j * i), ]
  #     grad = grad - t(as.matrix((p - y_ind)))%*%Xt ###negative of gradient due to computation of hessian below
  #   }
  #   return(grad)
  # }
  
  ##################################################################################################
  ###USE SIMPLIFIED TUNING HERE#####################################################################
  ##################################################################################################
  
  cat("initializing Metropolis candidate densities for ", nlgt, 
      " units ...", fill = TRUE)
  #Initialize betastar 
  betastarinit = c(rep(0, nvar))
  #Loop over individuals in the sample
  for (i in 1:nlgt) {
    if(length(lgtdata[[i]]$y)>=60){
      lgtdata[[i]] = c(lgtdata[[i]], list(converge = 0, 
                                          betafmle = c(rep(0, nvar)), hess = 10*diag(nvar)))
    }else{
      lgtdata[[i]] = c(lgtdata[[i]], list(converge = 0, 
                                          betafmle = c(rep(0, nvar)), hess = 1*diag(nvar)))
    }
    oldbetas[i, ] = lgtdata[[i]]$betafmle
    if (i%%50 == 0) 
      cat("  completed unit #", i, fill = TRUE)
    #fsh()
  }
  betas_tuning = oldbetas
  betastarpooled = betastarinit
  
  ###################
  ###END OF TUNING###
  ###################
  ind = NULL
  ninc = floor(nlgt/ncomp)
  for (i in 1:(ncomp - 1)) {
    ind = c(ind, rep(i, ninc))
  }
  if (ncomp != 1) {
    ###if there is more than a single component, equally divide individuals across components
    ind = c(ind, rep(ncomp, nlgt - length(ind)))
  }
  else {
    ind = rep(1, nlgt)
  }
  oldprob = rep(1/ncomp, ncomp)
  if (drawdelta) {
    olddelta = rep(0, nz * nvar)
  }
  else {
    olddelta = 0
    Z = matrix(0)
    deltabar = 0
    Ad = matrix(0)
  }
  draws = rhierMnlRwMixture_rcpp_loop_mod(lgtdata, nvar_c, Z, deltabar, Ad, mustarbarc, Amu, nu, V, Gammabar, A_Gamma, nu_Sigma, V_Sigma,
                                          s, R, keep, nprint, drawdelta, as.matrix(olddelta), a, oldprob, oldbetas, ind, betas_tuning, 
                                          betastarpooled, ind_stay_temp, ind_stay,ncomp)
  
  if (drawdelta) {
    attributes(draws$Deltadraw)$class = c("bayesm.mat", "mcmc")
    attributes(draws$Deltadraw)$mcpar = c(1, R, keep)
  }
  attributes(draws$betadraw)$class = c("bayesm.hcoef")
  attributes(draws$nmix)$class = "bayesm.nmix"
  return(draws)
}
