### Include compiled functions into rhierMnlRwMixture 
### NEW rhierMnlRwMixture function with varying choice sets and parts in CPP for optimal speed

rhierMnlRwMixture_varyCS_cpp <- function (Data, Prior, Mcmc) 
{
  pandterm=function(message) {stop(message,call.=FALSE)}
  
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
  
  
  if (missing(Data)) {
    pandterm("Requires Data argument -- list of p,lgtdata, and (possibly) Z")
  }
  # if (is.null(Data$p)) { # need p vector
  # pandterm("Requires Data element p (# chce alternatives)")
  # }
  # p = Data$p
  # 
  # if (is.null(Data$lgtdata)) {
  # pandterm("Requires Data element lgtdata (list of data for each unit)")
  # }
  lgtdata = Data$lgtdata
  nlgt = length(lgtdata)
  drawdelta = TRUE
  if (is.null(Data$Z)) {
    cat("Z not specified", fill = TRUE)
    fsh()
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
  }
  else {
    ncomp = Prior$ncomp
  }
  if (is.null(Prior$mubar)) {
    mubar = matrix(rep(0, nvar), nrow = 1)
  }
  else {
    mubar = matrix(Prior$mubar, nrow = 1)
  }
  if (ncol(mubar) != nvar) {
    pandterm(paste("mubar must have ncomp cols, ncol(mubar)= ", 
                   ncol(mubar)))
  }
  if (is.null(Prior$Amu)) {
    Amu = matrix(0.01, ncol = 1)
  }
  else {
    Amu = matrix(Prior$Amu, ncol = 1)
  }
  if (ncol(Amu) != 1 | nrow(Amu) != 1) {
    pandterm("Am must be a 1 x 1 array")
  }
  if (is.null(Prior$nu)) {
    nu = nvar + 5
    #nu = nvar + 15
    ###ALLENBY ET AL (2014)
  }
  else {
    nu = Prior$nu
  }
  if (nu < 1) {
    pandterm("invalid nu value")
  }
  #if(is.null(Prior$Ind_Restricted)){
  ####UPDATE
  #indrestr=6:28
  #Ind_Restricted = vector(mode="numeric",nvar)
  #Ind_Restricted[indrestr]=.5
  #Ind_Restricted[-indrestr]=1
  ###ALLENBY ET AL (2014)
  #}
  #else{
  #Ind_Restricted = Prior$Ind_Restricted
  #}
  if (is.null(Prior$V)) {
    #V = nu * Ind_Restricted * diag(nvar)
    V = nu * diag(nvar)
  }
  else {
    V = Prior$V
  }
  if (sum(dim(V) == c(nvar, nvar)) != 2) 
    pandterm("Invalid V in prior")
  if (is.null(Prior$Ad) & drawdelta) {
    Ad = 0.01 * diag(nvar * nz)
  }
  else {
    Ad = Prior$Ad
  }
  if (drawdelta) {
    if (ncol(Ad) != nvar * nz | nrow(Ad) != nvar * nz) {
      pandterm("Ad must be nvar*nz x nvar*nz")
    }
  }
  if (is.null(Prior$deltabar) & drawdelta) {
    deltabar = rep(0, nz * nvar)
  }
  else {
    deltabar = Prior$deltabar
  }
  if (drawdelta) {
    if (length(deltabar) != nz * nvar) {
      pandterm("deltabar must be of length nvar*nz")
    }
  }
  if (is.null(Prior$a)) {
    a = rep(5, ncomp)
  }
  else {
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
  }
  else {
    if (is.null(Mcmc$s)) {
      s = 2.38/sqrt(nvar)
    }
    else {
      s = Mcmc$s
    }
    if (is.null(Mcmc$w)) {
      w = 0.1
    }
    else {
      w = Mcmc$w
    }
    if (is.null(Mcmc$keep)) {
      keep = 1
    }
    else {
      keep = Mcmc$keep
    }
    if (is.null(Mcmc$R)) {
      pandterm("Requires R argument in Mcmc list")
    }
    else {
      R = Mcmc$R
    }
    if (is.null(Mcmc$nprint)) {
      nprint = 100
    }
    else {
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
  #    fill = TRUE)
  cat(paste("   for ", nlgt, " cross-sectional units"), fill = TRUE)
  cat(" ", fill = TRUE)
  cat("Prior Parms: ", fill = TRUE)
  cat("nu =", nu, fill = TRUE)
  cat("V ", fill = TRUE)
  print(V)
  cat("mubar ", fill = TRUE)
  print(mubar)
  cat("Amu ", fill = TRUE)
  print(Amu)
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
  
  #mnlHess_mod = function (betastar, y, X) {
  mnlHess_varyCS = function(beta,y,X,p) {
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
    n = length(y)
    k=ncol(X)
    Hess=matrix(double(k*k),ncol=k)
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
  
  ##################################################################################################
  ###TUNE CANDIDATES USING THE APPROXIMATIVE HESSIAN FOR BETASTARS USING DELTA METHOD###############
  ##################################################################################################
  cat("initializing Metropolis candidate densities for ", nlgt, 
      " units ...", fill = TRUE)
  fsh()
  
  betainit=c(rep(0,nvar))
  #
  #  compute pooled optimum
  #
  out=optim(par=betainit, fn=llmnl_varyCS,method="BFGS",control=list(fnscale=-1,trace=0,reltol=1e-6), 
            y=ypooled,X=Xpooled,choice=choicepooled,p=ppooled)
  # fnscale=-1 becaue we want to maximize
  # Used the Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm
  
  betapooled=out$par
  H=mnlHess_varyCS(beta=betapooled,y=ypooled,X=Xpooled,p=ppooled)
  rootH=chol(H)
  for (i in 1:nlgt) 
  {
    wgt=length(lgtdata[[i]]$y)/length(ypooled)
    out=optim(par=betapooled,fn=llmnlFract,method="BFGS",control=list( fnscale=-1,trace=0,reltol=1e-4), 
              X=lgtdata[[i]]$X,y=lgtdata[[i]]$y,choice=lgtdata[[i]]$choice,p=lgtdata[[i]]$p,betapooled=betapooled,rootH=rootH,w=w,wgt=wgt)
    if(out$convergence == 0) 
    { hess=mnlHess_varyCS(out$par,lgtdata[[i]]$y,lgtdata[[i]]$X,p=lgtdata[[i]]$p)
    lgtdata[[i]]=c(lgtdata[[i]],list(converge=1,betafmle=out$par,hess=hess)) }
    else
    { lgtdata[[i]]=c(lgtdata[[i]],list(converge=0,betafmle=c(rep(0,nvar)),
                                       hess=diag(nvar))) }
    oldbetas[i,]=lgtdata[[i]]$betafmle
    if(i%%50 ==0) cat("  completed unit #",i,fill=TRUE)
    #fsh()
  }
  
  betas_tuning = oldbetas
  ###################
  ###END OF TUNING###
  ###################
  ind = NULL
  ninc = floor(nlgt/ncomp)
  for (i in 1:(ncomp - 1)) {
    ind = c(ind, rep(i, ninc))
  }
  if (ncomp != 1) {
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
  #draws = rhierMnlRwMixture_rcpp_loop_mod(lgtdata, Z, deltabar, 
  #Ad, mubar, Amu, nu, V, s, R, keep, nprint, drawdelta, 
  #as.matrix(olddelta), a, oldprob, oldbetas, ind, betas_tuning, betastarpooled,
  #ind_stay_temp, ind_stay)
  draws = rhierMnlRwMixture_rcpp_loop_mod(lgtdata, Z, deltabar, 
                                          Ad, mubar, Amu, nu, V, s, R, keep, nprint, drawdelta, 
                                          as.matrix(olddelta), a, oldprob, oldbetas, ind, betas_tuning,
                                          ind_stay_temp, ind_stay)
  
  if (drawdelta) {
    attributes(draws$Deltadraw)$class = c("bayesm.mat", "mcmc")
    attributes(draws$Deltadraw)$mcpar = c(1, R, keep)
  }
  attributes(draws$betadraw)$class = c("bayesm.hcoef")
  attributes(draws$nmix)$class = "bayesm.nmix"
  return(draws)
}
