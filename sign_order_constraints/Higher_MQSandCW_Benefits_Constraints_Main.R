############################################################
# Replication files
# Benefits of constraints
# Marco Kotschedoff & Max Pachali
# March 2019
# Illustration: Benefits of prior constraints
# Full information choice set & seasonality controls

#################################################################################################################################################
### Note: This replication code is based on a randomly drawn subsample of N = 119 households to reduce computation time #########################
### Results presented in slides are based a 10x bigger sample of N = 1337 households and therefore not exactly recovered in this illustration ###      
#################################################################################################################################################

############################################################
# Set up
############################################################

#remove everything
rm(list=ls())

library(foreign)
library(bayesm)
library(Rcpp)
library(devtools)
library(RcppArmadillo)
library(matrixStats)
library(Matrix)
library(MASS)
library(ggplot2)   
library(plyr)
library(stargazer)
library(doBy)
library(pracma)
library(xtable)
library(tikzDevice)
library(plyr)
library(latex2exp)
library(coda)

### Increase memory capacities
memory.limit(size=100000000)

# ### Path settings - starting from path of replications files
getwd()
setwd("C:\\Users\\pachali\\PhD_Economics\\Sign_And_Order_Case_Illus\\Egg_Case_Replication")
setwd("YOUR PATH TO REPLICATION FILES")

### Load replication data & demographics
load("Synthetic_Replication_Data_1000.Rdata")

### Inspect data
N = length(lgtdata_sub)
#View(lgtdata_sub[[1]]$X)
#lgtdata_sub[[1]]$y
#lgtdata_sub[[1]]$p

########################################################################################################
###Prepare estimation data for MCMC inference: Re-ordering necessary for constrained sampler later on###
########################################################################################################

### Permutate the Design matrix to separate constrained & unconstrained coefficients: We only use this to conveniently compare posterior outcomes later on...
nvar = dim(lgtdata_sub[[1]]$X)[2] - 1 ###(!!!identification restriction!!!)
n_stores = 10 ###number of stores included
nvar_uc = 3 + n_stores - 1 #constrain coefficients except for cage as well as store parameters (!!!identification restriction!!!) 
nvar_c = nvar-nvar_uc 

ind_uc = c(0,1,rep(0,4),rep(1,n_stores),1,1,rep(0,9)) #equals 1 if a coefficient is UNconstrained and 0 otherwise
seqdimX = seq(dim(lgtdata_sub[[1]]$X)[2])
lgtdata_c = NULL

for(i in 1:N){
  X_uc = lgtdata_sub[[i]]$X[,ind_uc*seqdimX]
  X_c = lgtdata_sub[[i]]$X[,-ind_uc*seqdimX]
  X_temp = cbind(X_c,X_uc)
  ident = which(colnames(X_temp)=="RT 1")
  lgtdata_c[[i]] = list(y = lgtdata_sub[[i]]$y, X = X_temp[,-ident]
                        , p = lgtdata_sub[[i]]$p, hhid = lgtdata_sub[[i]]$hhid, choice = lgtdata_sub[[i]]$choice, nobs = lgtdata_sub[[i]]$nobs) ###set identification restriction here: first store is baseline
}

#View(lgtdata_c[[1]]$X)

###################################################################
### Specify standard priors for NM1-model as in Rossi et al. (2005)
###################################################################

# Unique components
ncomp = NULL
ncomp[[1]] = 1

# Mcmc -depend whether these are the smaller mock data
keep = 100
R = 400000
mcmc1 = list(keep=keep,R=R,w=0.1) 

#########################################################
### Run MCMC ############################################
#########################################################
### Compile MCMC sampler
Rcpp::sourceCpp('rhierMnlRwMixture_rcpp_loop_varyCS_unconstrained.cpp',showOutput = FALSE)
source('rhierMnlRwMixture_varyCS_unconstrained_rcpp.R')

out_UC = rhierMnlRwMixture_varyCS_cpp(Data=list(lgtdata=lgtdata_c),Prior=list(ncomp = ncomp[[1]]), Mcmc=mcmc1)

#########################################
###Specify about burn-in#################
#########################################
R = dim(out_UC$betadraw)[3]

burnin = 2000 

betastar_UC = out_UC$betadraw
compdraw_UC = out_UC$nmix$compdraw
probdraw_UC = out_UC$nmix$probdraw
rejection_UC = out_UC$rejection
loglike_UC=out_UC$loglike

windows()
plot(loglike_UC, type = "l")

### Get rid of burnin
betastar_UC = betastar_UC[,,(burnin+1):R]
compdraw_UC = compdraw_UC[(burnin+1):R]
probdraw_UC = probdraw_UC[(burnin+1):R]
rejection_UC = rejection_UC[(burnin+1):R,]
loglike_UC = loglike_UC[(burnin+1):R]

R = dim(betastar_UC)[3]
N = dim(betastar_UC)[1]
  
#############################
###Log marginal likelihood###
#############################
one_perc = R*0.01
loglike_trimmed_UC = sort(loglike_UC)[(one_perc+1):(R-one_perc-1)]
# NRE_trimmed:
NRE_trimmed_UC = logMargDenNR(loglike_trimmed_UC)


#######################################################
###Analyze Posterior Preference Distribution###########
#######################################################
# !!!Reobtain previous ordering: do this only once!!!
colnames(betastar_UC) = colnames(lgtdata_c[[1]]$X)
betastar_UC = betastar_UC[,c(1,15,2,3,4,5,25,6,8,10,26,7,9,11,12,13,14,16:24),]
beta_UC = betastar_UC 

# Preference distribution based on llmns
set.seed(66)
l = 30
beta_sim_UC <- array(0,dim=c(R*l,dim(betastar_UC)[2]))
index=rep(rank(runif(R)),l)
index_n = rep(rank(runif(N)),round((R*l)/N)+1)

# Data generation
for(j in 1:(R*l)){
  beta_sim_UC[j,] = beta_UC[index_n[j],,index[j]]
}

colnames(beta_sim_UC) = colnames(beta_UC)

summary(beta_sim_UC)


#########################################################################################
######RUN CONSTRAINED MCMC SAMPLER NOW###################################################
#########################################################################################

### Specify subjective priors for NM1-model as in Appendix A.3
mustarbarc = matrix(rep(0, nvar_c), nrow = nvar_c)
Gammabar = matrix(0, nrow = (nvar_c+1), ncol = nvar_uc)
A_Gamma = diag((1/100), nrow = (nvar_c+1), ncol = (nvar_c+1))
nu_Sigma = nvar_uc + 5 
V_Sigma = nu_Sigma * diag(nvar_uc)

# Amu of constrained coefficients 
Amustar_entry = diag(nvar_c)*1
#Price
Amustar_entry[1,1] = 1/4 #regular
#Six
Amustar_entry[5,5] = 1/4
#Barn 
Amustar_entry[2,2] = 1/2 #regular
Amustar_entry[6,6] = 1/2 #Eastern
Amustar_entry[7,7] = 1/2 #Xmas
Amustar_entry[12,12] = 1/2 #post ban
#Free-Range
Amustar_entry[3,3] = 1/2 #regular
Amustar_entry[8,8] = 1/2 #Eastern
Amustar_entry[9,9] = 1/2 #Xmas
Amustar_entry[13,13] = 1/2 #post ban
#Organic
Amustar_entry[4,4] = 1 #regular
Amustar_entry[10,10] = 1 #Eastern
Amustar_entry[11,11] = 1 #Xmas
Amustar_entry[14,14] = 1 #post ban

Amu = NULL 
Amu[[1]] = Amustar_entry

nu = NULL
nu[[1]] = 30

V = NULL
V[[1]] = nu[[1]] * diag(nvar_c)*0.25

Z = NULL
Ad = 0.01 * diag(nvar * ncol(Z))
ncomp = NULL
ncomp[[1]] = 1

Model_settings = list(nvar=nvar,Amu=Amu,Ad=Ad,V=V,ncomp=ncomp, nu=nu)
keep = 100
R = 400000
mcmc1 = list(keep=keep,R=R,w=0.1,s=c(rep(0.5,nvar_c),rep(0.8,nvar_uc))) # vary step size if needed

### Run MCMC 
# Compile MCMC sampler
Rcpp::sourceCpp('rhierMnlRwMixture_rcpp_loop_eggs_MCD_untuned_controls_full.cpp',showOutput = FALSE)
source('rhierMnlRwMixture_main_eggs_untuned.R')

scen = 1
out_C = rhierMnlRwMixture_SR(Data=list(lgtdata=lgtdata_c),Prior=list(ncomp = ncomp[[scen]], nu = nu[[scen]], V = V[[scen]], Amu = Amu[[scen]],
                                                                      mustarbarc = mustarbarc, Gammabar = Gammabar, A_Gamma = A_Gamma, nu_Sigma = nu_Sigma,
                                                                      V_Sigma = V_Sigma),Mcmc=mcmc1, nvar_c = nvar_c)

#####################
### Specify burnin###
#####################
R = dim(out_C$betadraw)[3]
burnin = 2000 

betastar_C = out_C$betadraw
compdraw_C = out_C$nmix$compdraw
probdraw_C = out_C$nmix$probdraw
rejection_C = out_C$rejection
loglike_C = out_C$loglike

windows()
plot(loglike_C, type="l")

betastar_C = betastar_C[,,(burnin+1):R]
compdraw_C = compdraw_C[(burnin+1):R]
if(scen==1){
  probdraw_C = probdraw_C[(burnin+1):R]
}else{
  probdraw_C = probdraw_C[(burnin+1):R,]
}
rejection_C = rejection_C[(burnin+1):R,]
loglike_C = loglike_C[(burnin+1):R]

# Update R
R = dim(betastar_C)[3]
N = dim(betastar_C)[1]

##################################################################################################
### Compute Log marginal likelihood as in Table 7 (note that results differ due to much smaller N)
##################################################################################################
one_perc = R*0.01
loglike_trimmed_C = sort(loglike_C)[(one_perc+1):(R-one_perc-1)]
NRE_trimmed_C = logMargDenNR(loglike_trimmed_C)
# Compare Constrained & Unconstrained
NRE_trimmed_C
NRE_trimmed_UC

#######################################################################
### Derive posterior distribution of household preferences#############
#######################################################################

### Transform betastars to betas & reobtain previous ordering
colnames(betastar_C) = colnames(lgtdata_c[[1]]$X)
betastar_C = betastar_C[,c(1,15,2,3,4,5,25,6,8,10,26,7,9,11,12,13,14,16:24),]
beta_C = betastar_C 
# Price 
beta_C[,1,] = -exp(betastar_C[,1,])
# Cage up to organic eggs
beta_C[,2,] = beta_C[,2,] #unconstrained
beta_C[,3,] = beta_C[,2,] + exp(betastar_C[,3,]) #barn
beta_C[,4,] = beta_C[,3,] + exp(betastar_C[,4,]) #free range
beta_C[,5,] = beta_C[,4,] + exp(betastar_C[,5,]) #organic
# Six relative to ten
beta_C[,6,] = -exp(betastar_C[,6,])
# Cage up to organic eggs in Eastern
beta_C[,7,] = beta_C[,7,] #unconstrained
beta_C[,8,] = beta_C[,7,] + exp(betastar_C[,8,]) #barn
beta_C[,9,] = beta_C[,8,] + exp(betastar_C[,9,]) #free range
beta_C[,10,] = beta_C[,9,] + exp(betastar_C[,10,]) #organic
# Cage up to organic eggs in Xmas
beta_C[,11,] = beta_C[,11,] #unconstrained
beta_C[,12,] = beta_C[,11,] + exp(betastar_C[,12,]) #barn
beta_C[,13,] = beta_C[,12,] + exp(betastar_C[,13,]) #free range
beta_C[,14,] = beta_C[,13,] + exp(betastar_C[,14,]) #organic
# Barn up to organic eggs after regime change
beta_C[,15,] = beta_C[,2,] + exp(betastar_C[,15,]) #barn
beta_C[,16,] = beta_C[,15,] + exp(betastar_C[,16,]) #free range
beta_C[,17,] = beta_C[,16,] + exp(betastar_C[,17,]) #organic

# create preference distribution
beta_sim_C <- array(0,dim=c(R*l,dim(betastar_C)[2]))
# Data generation
for(j in 1:(R*l)){
  beta_sim_C[j,] = beta_C[index_n[j],,index[j]]
}

colnames(beta_sim_C) = colnames(beta_C)

### Compare the two:
# Unconstrained
summary(beta_sim_UC[,1:6])
# Constrained
summary(beta_sim_C[,1:6])

#####################################################################################
### Compare plots of marginal posteriors b/w constrained & unconstrained version#####
#####################################################################################
# price coefficient
pricecoeff_frame = data.frame(Preferences = as.vector(c(beta_sim_UC[,1],beta_sim_C[,1])), Form = c(rep("Unconstrained",length(beta_sim_UC[,1])),
                                                                                                  rep("Constrained",length(beta_sim_C[,1]))))

pricecoeff_frame$Form <- factor(pricecoeff_frame$Form,levels=c("Unconstrained","Constrained"))

windows()
ggplot(pricecoeff_frame, aes(Preferences, fill = Form, colour = Form)) +
  geom_density(alpha = 0.1, size=1.2, bw=0.3, aes(linetype = Form))+
  #geom_line(size=1.5,aes(linetype = model)) +
  xlim(-15,5) +
  ylim(0,0.25) +
  xlab("coefficient") +
  ylab("density") +
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),  
        axis.title.x = element_text(size=22),
        #axis.title.y = element_text(size=22),legend.title = element_blank(),legend.text=element_blank(),legend.position="none") +
        axis.title.y = element_text(size=22),legend.title = element_blank(),legend.text=element_text(size=16),legend.justification = c(0, 1),
        legend.position = c(0, 1), legend.background = element_rect(fill = "transparent", linetype = "dotted")) +
  scale_fill_manual(values = c("brown4","brown4")) +
  scale_color_manual(values = c("brown4","brown4"))  

# Battery Uncon vs. Constrained
coeff_frame = data.frame(Preferences = as.vector(c(beta_sim_UC[,2],beta_sim_C[,2])), Form = c(rep("Unconstrained",length(beta_sim_UC[,2])),
                                                rep("Constrained",length(beta_sim_C[,2]))))

coeff_frame$Form <- factor(coeff_frame$Form,levels=c("Unconstrained","Constrained"))

windows()
ggplot(coeff_frame, aes(Preferences, fill = Form, colour = Form)) +
  geom_density(alpha = 0.1,size=1.2,bw=2.5, aes(linetype = Form))+
  xlim(-25,20) +
  ylim(0,0.1) +
  xlab("coefficient") +
  ylab("density") +
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),  
        axis.title.x = element_text(size=22),
        axis.title.y = element_text(size=22),legend.title = element_blank(),legend.text=element_text(size=16),legend.justification = c(0, 1),
        legend.position = c(0, 1), legend.background = element_rect(fill = "transparent", linetype = "dotted")) +
  scale_fill_manual(values = c("blue","blue")) +
  scale_color_manual(values = c("blue","blue"))  

# Organic Uncon vs. Constrained
coeff_frame = data.frame(Preferences = as.vector(c(beta_sim_UC[,5],beta_sim_C[,5])), Form = c(rep("Unconstrained",length(beta_sim_UC[,5])),
                                                                                              rep("Constrained",length(beta_sim_C[,5]))))

coeff_frame$Form <- factor(coeff_frame$Form,levels=c("Unconstrained","Constrained"))

windows()
ggplot(coeff_frame, aes(Preferences, fill = Form, colour = Form)) +
  geom_density(alpha = 0.1,size=1.2,bw=2.5, aes(linetype = Form))+
  xlim(-25,20) +
  ylim(0,0.1) +
  xlab("coefficient") +
  ylab("density") +
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),  
        axis.title.x = element_text(size=22),
        axis.title.y = element_text(size=22),legend.title = element_blank(),legend.text=element_text(size=16),legend.justification = c(0, 1),
        legend.position = c(0, 1), legend.background = element_rect(fill = "transparent", linetype = "dotted")) +
  scale_fill_manual(values = c("black","black")) +
  scale_color_manual(values = c("black","black"))  

###################################################
###########Impute Marginal Costs###################
###################################################
load("X_markets_etc.Rdata")
X_market=X_markets[,-7] ###Identification restriction

### Source counterfactual functions
source('Functions_Counterfactual_withstoresMK_rccp.R')
Rcpp::sourceCpp('Speed++_Logit_probas.cpp',showOutput = FALSE)

### Some settings for numerical analysis
maxiter_set=200000
tol_set=0.1^6

### Keep only one outside good
out_ind=which(X_market$product=="Outside option 10 ST", arr.ind = TRUE)
X_market[out_ind,]$retailer_anonym=c("Outside")

### Ownerschip matrix real: real ownership: e.g. penny and rewe belong to rewe-group
ownership_real=diag(dim(X_market)[1])
for (col in 1:dim(X_market)[1]) {
  for (row in 1:dim(X_market)[1]) {
    if (X_market$chain_anonym[col]==X_market$chain_anonym[row]){ownership_real[row,col]=1}
  }
}

# Outside good: put everywhere zeros!
out_ind=which(X_market$product=="Outside option 10 ST", arr.ind = TRUE)
ownership_real[out_ind,]=0
ownership_real[,out_ind]=0

# Update betas
beta_sim_UC=cbind(beta_sim_UC[,1:6],beta_sim_UC[,18:26]) #counter-factuals w/o regime change & off-season
beta_sim_C=cbind(beta_sim_C[,1:6],beta_sim_C[,18:26]) #counter-factuals w/o regime change & off-season

# Update market shares given preference, because now only 10 ST products
ms_real_UC=logit_probas(beta_draws=beta_sim_UC,X=X_market,prod_names=X_market$product)$mean
ms_real_C=logit_probas(beta_draws=beta_sim_C,X=X_market,prod_names=X_market$product)$mean

X_market_UC = X_market; X_market_UC$ms <- ms_real_UC 
X_market_C = X_market; X_market_C$ms <- ms_real_C 

#Xdata_UC=X_market_UC; Xdata_C=X_market_C; 

MC_real_UC=MC_logit_NE(beta_draws=beta_sim_UC,X=X_market_UC,price_vec=X_market_UC$price,position=1,ownership=ownership_real,prod_names=X_market_UC$product,
                    firm_names=X_market_UC$retailer_anonym)

MC_real_C=MC_logit_NE(beta_draws=beta_sim_C,X=X_market_C,price_vec=X_market_C$price,position=1,ownership=ownership_real,prod_names=X_market_C$product,
                       firm_names=X_market_C$retailer_anonym)


MC_real_full_UC = MC_real_UC 
MC_real_full_C = MC_real_C 


firm_costs_UC=data.frame(MC=MC_real_UC$MC,price=X_market_UC$price[2:length(X_market_UC$price)],ms=X_market_UC$ms[2:length(X_market_UC$ms)],
                         product=MC_real_UC$products,firm=MC_real_UC$firms,margin=MC_real_UC$Margin)

firm_costs_C=data.frame(MC=MC_real_C$MC,price=X_market_C$price[2:length(X_market_C$price)],ms=X_market_C$ms[2:length(X_market_C$ms)],
                         product=MC_real_C$products,firm=MC_real_C$firms,margin=MC_real_C$Margin)


### Compare market-share weighted marginal cost estimates of constrained & unconstrained models
weighted_mean_cost_UC = ddply(firm_costs_UC,.(product),summarise, MC = weighted.mean(MC,ms))

weighted_mean_cost_C = ddply(firm_costs_C,.(product),summarise, MC = weighted.mean(MC,ms))

Order_costs = c(4,2,3,1)
Compare_costs = cbind(weighted_mean_cost_UC[Order_costs,2],weighted_mean_cost_C[Order_costs,2])

rownames(Compare_costs) <- c("Battery 10 units","Barn 10 units","Free-range 10 units", "Organic 10 units")
colnames(Compare_costs) <- c("Unconstrained","Constrained")

#################################################################################################################################################
### Note: This replication code is based on a randomly drawn subsample of N = 119 households to reduce computation time #########################
### Results presented in slides are based a 10x bigger sample of N = 1337 households and therefore not exactly recovered in this illustration ###      
#################################################################################################################################################
Compare_costs





