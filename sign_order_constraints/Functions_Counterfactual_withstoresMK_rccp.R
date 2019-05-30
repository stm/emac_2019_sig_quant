#######################################################
# A collection of functions for Counterfactual Scenarios
# - Market shares
# - Profits
# - FOCs
# - Compensation Variation - Change in Consumer Welfare
# etc.

# Written by Marco Kotschedoff & Max Pachali Jan 2017
######################################################




logit_CV = function(beta_draws,X,X_counter,market_size) {
  
  #This is the stabilized version of the computing changes 
  # in Consumer welfare.
  
  #Define X-matrix
  X=as.matrix(X[,1:dim(beta_draws)[2]])
  X_counter=as.matrix(X_counter[,1:dim(beta_draws)[2]])
  
  #Compute Indirect utilities
  V=X%*%t(beta_draws)
  V_counter=X_counter%*%t(beta_draws)
  
  #Stabilize by substracting maximum V (of both scenarios) of each draw (does not change probas)
  # It is thus a stabilizing constant
  V.col.max <- apply(V, 2, max)
  V_counter.col.max <- apply(V_counter, 2, max)
  V.col.max_combine<- apply(rbind(V.col.max,V_counter.col.max),2,max)
  #now substract vector from this matrix but from each element by column
  V_stab=sweep(V,MARGIN=2,V.col.max_combine,FUN="-")
  V_counter_stab=sweep(V_counter,MARGIN=2,V.col.max_combine,FUN="-")
  
  #Compute Exponential
  e_Xb=t(exp(V_stab))
  e_Xb_counter=t(exp(V_counter_stab))
  
  #Compute log sums
  log_sum=log(rowSums(e_Xb))
  log_sum_counter=log(rowSums(e_Xb_counter))
  
  
  #Compute compensation variation
  CV_draws=1/abs(beta_draws[,1])*(log_sum_counter-log_sum)
  
  #Compute compensation variation if only determinsitic part counts
  CV_draws_deterministic=1/abs(beta_draws[,1])*(V_counter.col.max-V.col.max)
  
  
  #Alternative but identical #Delta_CS=market_size*mean(CV_draws)
  Delta_CS=sum(market_size/length(CV_draws)*CV_draws)
  Delta_CS_deterministic=sum(market_size/length(CV_draws)*CV_draws_deterministic)
  pos_ind=which(CV_draws>0, arr.ind = TRUE)
  neg_ind=which(CV_draws<0, arr.ind = TRUE)
  pos_ind_deterministic=which(CV_draws_deterministic>0, arr.ind = TRUE)
  neg_ind_deterministic=which(CV_draws_deterministic<0, arr.ind = TRUE)
  Delta_CS_losers=sum((market_size/length(CV_draws))*CV_draws[neg_ind])
  Delta_CS_winners=sum((market_size/length(CV_draws))*CV_draws[pos_ind])
  Delta_CS_losers_deterministic=sum((market_size/length(CV_draws))*CV_draws_deterministic[neg_ind_deterministic])
  Delta_CS_winners_deterministic=sum((market_size/length(CV_draws))*CV_draws_deterministic[pos_ind_deterministic])
  #Results:
  logit_CV<-NULL
  logit_CV$CV_draws<-CV_draws
  logit_CV$Delta_CS<-Delta_CS
  logit_CV$Delta_CS_losers<-Delta_CS_losers
  logit_CV$Delta_CS_winners<-Delta_CS_winners
  #Deterministic
  logit_CV$CV_draws_deterministic<-CV_draws_deterministic
  logit_CV$Delta_CS_deterministic<-Delta_CS_deterministic
  logit_CV$Delta_CS_losers_deterministic<-Delta_CS_losers_deterministic
  logit_CV$Delta_CS_winners_deterministic<-Delta_CS_winners_deterministic
  return(logit_CV)
}


logit_CV_explicit = function(beta_draws,X,X_counter,market_size, extreme_value_draws, counter_ind) {
  
  #This is the stabilized version of the computing changes 
  # in Consumer welfare (explicit computation)
  # 
  
  #Define X-matrix
  X=as.matrix(X[,1:dim(beta_draws)[2]])
  X_counter=as.matrix(X_counter[,1:dim(beta_draws)[2]])
  
  #Compute Indirect utilities: only deterministic part
  V=X%*%t(beta_draws)
  V_counter=X_counter%*%t(beta_draws)
  
  #Add extreme value type one draws here explicity: Get indirect utility
  U=V+extreme_value_draws
  U_counter=V_counter+extreme_value_draws[counter_ind,]
  
  #Do not need to stabilze because we do use the log or exponential function
  # Calc. maxiumum utlitity
  U.col.max <- as.vector(apply(U, 2, max))
  U_counter.col.max <-  as.vector(apply(U_counter, 2, max))

  #Compute compensation variation with explicit extreme value type I draws
  CV_draws=1/abs(beta_draws[,1])*(U_counter.col.max-U.col.max)

  
  #Alternative but identical #Delta_CS=market_size*mean(CV_draws)
  Delta_CS=sum(market_size/length(CV_draws)*CV_draws)
  pos_ind=which(CV_draws>0, arr.ind = TRUE)
  neg_ind=which(CV_draws<0, arr.ind = TRUE)
  Delta_CS_losers=sum((market_size/length(CV_draws))*CV_draws[neg_ind])
  Delta_CS_winners=sum((market_size/length(CV_draws))*CV_draws[pos_ind])
  #Results:
  logit_CV<-NULL
  logit_CV$CV_draws<-CV_draws
  logit_CV$Delta_CS<-Delta_CS
  logit_CV$Delta_CS_losers<-Delta_CS_losers
  logit_CV$Delta_CS_winners<-Delta_CS_winners
  return(logit_CV)
}


logit_CV_old = function(beta_draws,X,X_counter,market_size) {
  
  #Old version - not stabilized
  
  #status quo
  X=as.matrix(X[,1:dim(beta_draws)[2]])
  e_Xb=t(exp(X%*%t(beta_draws)))
  log_sum=log(rowSums(e_Xb))
  
  #counterfactual
  X_counter=as.matrix(X_counter[,1:dim(beta_draws)[2]])
  e_Xb_counter=t(exp(X_counter%*%t(beta_draws)))
  log_sum_counter=log(rowSums(e_Xb_counter))
  
  CV_draws=1/abs(beta_draws[,1])*(log_sum_counter-log_sum)
  
  
  #Alternative but identical #Delta_CS=market_size*mean(CV_draws)
  Delta_CS=sum(market_size/length(CV_draws)*CV_draws)
  pos_ind=which(CV_draws>0, arr.ind = TRUE)
  neg_ind=which(CV_draws<0, arr.ind = TRUE)
  Delta_CS_losers=sum((market_size/length(CV_draws))*CV_draws[neg_ind])
  Delta_CS_winners=sum((market_size/length(CV_draws))*CV_draws[pos_ind])
  #Results:
  logit_CV<-NULL
  logit_CV$CV_draws<-CV_draws
  logit_CV$Delta_CS<-Delta_CS
  logit_CV$Delta_CS_losers<-Delta_CS_losers
  logit_CV$Delta_CS_winners<-Delta_CS_winners
  return(logit_CV)
}


logit_profit= function(beta_draws,X,price_vec,cost_vec,market_size,prod_names,firm_names) {
  
  #replace prices in X with given price vector
  X$price=price_vec
  market_shares=logit_probas(beta_draws,X,prod_names)$mean
  
  #strip outside good from matrixes
  out_ind=which(X$product=="Outside option 10 ST", arr.ind = TRUE)
  market_shares=market_shares[-out_ind]
  price_vec=price_vec[-out_ind]
  profit=(price_vec-cost_vec)*market_shares*market_size
  results <- data.frame(profit=profit,firm_names=firm_names[-out_ind])
  #sum profits by firm
  firm_profits=aggregate(results$profit, by=list(Category=results$firm_names), FUN=sum)
  #Results:
  logit_profit<-NULL
  logit_profit$firm<-levels(as.factor(firm_names))
  logit_profit$profit<-firm_profits$x
  return(logit_profit)
}   






ban_FOC= function(price_vec,MC,ownership,Xdata,beta_draws,ban_vec) {
  
  #Define constant values here
  ban_ind=which(Xdata$product %in% ban_vec, arr.ind = TRUE)
  #second coefficient is battery egg
  
  X=Xdata[-ban_ind,]
  ownership=ownership[-ban_ind,-ban_ind]
  prod_names=X$product
  firm_names=X$retailer_anonym
  
  ban_ind_MC=which(MC$product %in% ban_vec, arr.ind = TRUE)
  cost_vec=MC$MC[-ban_ind_MC]
  position=1
  
  
  #FOC_result= logit_FOC_NE(beta_draws,X,price_vec,cost_vec,position,ownership,prod_names)
  FOC_result= logit_FOC_NE_markup(beta_draws,X,price_vec,cost_vec,position,ownership,prod_names)
  
  
  
  #Results:
  FOC_result
}





#find optimal subsidy
opt_subs= function(subsidy,MC,ownership,Xdata,beta_draws,min_stand,ban_vec) {
  
  min_stand_ind1=which(MC$product==min_stand[1], arr.ind = TRUE)
  min_stand_ind2=which(MC$product==min_stand[2], arr.ind = TRUE)
  MC_new=MC
  MC_new$MC[min_stand_ind1]=MC_new$MC[min_stand_ind1]-rep(subsidy*10, length(MC_new$MC[min_stand_ind1])) #since we consider only ten egg package sizes
  MC_new$MC[min_stand_ind2]=MC_new$MC[min_stand_ind2]-rep(subsidy*6, length(MC_new$MC[min_stand_ind2])) #six pack
  
  ban_ind=which(Xdata$product %in% ban_vec, arr.ind = TRUE)
  X_ban=Xdata[-ban_ind,]
  p0=X_ban$price
  p_new=fsolve(ban_FOC,p0,maxiter =maxiter_set, tol =tol_set,MC=MC_new,ownership=ownership,Xdata=Xdata,beta_draws=beta_draws,ban_vec=ban_vec)
  p_new=p_new$x
  X_ban$price<-p_new
  
  battery10_ind=which(Xdata$product=="KAEFIG 10 ST", arr.ind = TRUE)
  min_stand_ind1=which(X_ban$product==min_stand[1], arr.ind = TRUE)
  ms_preban=logit_probas(beta_draws,Xdata,prod_names=Xdata$product)$mean
  ms_postban=logit_probas(beta_draws,X_ban,prod_names=X_ban$product)$mean 
  Xdata$ms<-ms_preban
  X_ban$ms<-ms_postban
  # If there are different prices for the same product: 
  # use market share weighted mean as target
  weighted.mean.prices_preban=ddply(Xdata,.(product),summarise, price = weighted.mean(price,ms))
  weighted.mean.prices_postban=ddply(X_ban,.(product),summarise, price = weighted.mean(price,ms))
  price_min_stand_new=weighted.mean.prices_postban$price[weighted.mean.prices_postban$product==min_stand[1]]
  price_battery10_old=weighted.mean.prices_preban$price[weighted.mean.prices_preban$product=="KAEFIG 10 ST"]
  
  #Results:
  opt_subs=price_min_stand_new-price_battery10_old
  return(opt_subs)
}




#
#
logit_probas = function(beta_draws,X,prod_names) {
  
  #Compute Logit probas - stabilize by substracting maximum value
  # of indirect utility
  
  ###Transform X-matrix
  X = as.matrix(X[,1:dim(beta_draws)[2]])
  ###Initialize draws of market shares
  probas_draws=array(0,dim=c(dim(beta_draws)[1],dim(X)[1]))
  ###Compute stabilized market shares using rcpp code
  probas_draws = probabilities_Logit_cpp(beta_draws,X)
  
  probas.mean <- apply(probas_draws, c(2),mean)
  probas.median <- apply(probas_draws, c(2),median)
  
  #Results:
  logit_probas<-NULL
  logit_probas$draws=probas_draws
  logit_probas$mean=probas.mean
  logit_probas$median=probas.median
  logit_probas$products=prod_names
  return(logit_probas)
}






###    
logit_demand= function(beta_draws,X,market_size,prod_names) {
  
  
  market_shares=logit_probas(beta_draws,X,prod_names)$mean
  quantity=market_shares*market_size
  
  #Results:
  logit_demand<-NULL
  logit_demand$quantity=quantity
  logit_demand$market_shares=market_shares
  logit_demand$products=prod_names
  return(logit_demand)
}    


###    
logit_jacobian = function(beta_draws, X, position,prod_names) {
  if (is.null(position)) {
    position = 1
  }
  n_sample = dim(beta_draws)[1]
  n_products = dim(X)[1]
  logit_probas.draws = logit_probas(beta_draws,X,prod_names)$draws
  
  #Jacobian matrix for each individual
  jacobian_draws = array(0, dim = c(n_products, n_products, n_sample))
  
  for(draw in 1:n_sample){
    logit_probas_vec = as.vector(logit_probas.draws[draw, ])
    jacobian_draws[,, draw] = (-1) * beta_draws[draw, position] * tcrossprod(logit_probas_vec)
    diag(jacobian_draws[,, draw]) = beta_draws[draw, position] * logit_probas_vec * (1 - logit_probas_vec)
  }  
  
  # for (draw in 1:n_sample) {
  #   for (prd_j in 1:n_products) {
  #     jacobian_draws[prd_j, prd_j, draw] = beta_draws[draw, position] * logit_probas.draws[draw, prd_j] *
  #       (1 - logit_probas.draws[draw, prd_j])
  #     for (prd_k in 1:n_products) {
  #       if (prd_k != prd_j) {
  #         jacobian_draws[prd_j, prd_k, draw] = (-1) * beta_draws[draw, position] *
  #           logit_probas.draws[draw, prd_j] * logit_probas.draws[draw, prd_k]
  #       }
  #     }
  #   }
  # }
  
  jacobian.mean <- apply(jacobian_draws, c(1, 2), mean)
  jacobian.median <- apply(jacobian_draws, c(1, 2), median)
  
  #Results:
  logit_jacobian <- NULL
  logit_jacobian$draws = jacobian_draws
  logit_jacobian$mean = jacobian.mean
  logit_jacobian$median = jacobian.median
  logit_jacobian$products=prod_names
  return(logit_jacobian)
}  

###
logit_elasticity = function(beta_draws, X, position,prod_names) {
  if (is.null(position)) {
    position = 1
  }
  n_sample = dim(beta_draws)[1]
  n_products = dim(X)[1]
  logit_probas.draws = logit_probas(beta_draws,X,prod_names)$draws
  logit_probas.mean = logit_probas(beta_draws,X,prod_names)$mean
  logit_probas.median = logit_probas(beta_draws,X,prod_names)$median
  logit_jacobian.mean = logit_jacobian(beta_draws,X,position,prod_names)$mean
  logit_jacobian.median = logit_jacobian(beta_draws,X,position,prod_names)$median
  #Elasticity matrix for draws
  #Formular for mixed logit models
  #draws only for the jacobian necessary and then the integral-> can take the mean jacobiab
  elasticity = array(NaN, dim = c(n_products, n_products))
  
  #Mean elasticity
  for (prd_j in 1:n_products) {
    elasticity[prd_j, prd_j] = X[prd_j,position]/logit_probas.mean[prd_j]*logit_jacobian.mean[prd_j,prd_j]
    for (prd_k in 1:n_products) {
      if (prd_k != prd_j) {
        elasticity[prd_j, prd_k] = X[prd_k,position]/logit_probas.mean[prd_j]*logit_jacobian.mean[prd_k,prd_j]
      }
    }
  }
  #Median elasticity
  elasticity.median = array(NaN, dim = c(n_products, n_products))
  
  for (prd_j in 1:n_products) {
    elasticity.median[prd_j, prd_j] = X[prd_j,position]/logit_probas.median[prd_j]*logit_jacobian.median[prd_j,prd_j]
    for (prd_k in 1:n_products) {
      if (prd_k != prd_j) {
        elasticity.median[prd_j, prd_k] = X[prd_k,position]/logit_probas.median[prd_j]*logit_jacobian.median[prd_k,prd_j]
      }
    }
  }
  
  
  #Results:
  logit_elasiticty <- NULL
  logit_elasiticty$mean = elasticity
  logit_elasiticty$median = elasticity.median
  logit_elasiticty$products=prod_names
  return(logit_elasiticty)
}  

###


logit_FOC_NE= function(beta_draws,X,price_vec,cost_vec,position,ownership,prod_names) {
  
  #Indicate inside goods
  inside_row=which(rowSums(ownership) != 0, arr.ind = TRUE)
  #update prices, only inside goods
  X$price[inside_row]=price_vec[inside_row]
  
  
  #Compute probas given new price, other probas and jacobian do not change
  logit_probas.mean = logit_probas(beta_draws,X,prod_names)$mean
  logit_jacobian.mean = logit_jacobian(beta_draws,X,position,prod_names)$mean
  
  #strip outside good from matrixes
  strip_col=which(colSums(ownership) == 0, arr.ind = TRUE)
  strip_row=which(rowSums(ownership) == 0, arr.ind = TRUE)
  logit_jacobian.mean= logit_jacobian.mean[-strip_row,-strip_col]
  ownership_strip= ownership[-strip_row,-strip_col]
  logit_probas.mean=  logit_probas.mean[-strip_row]
  price_vec=price_vec[-strip_row]
  cost_vec=cost_vec
  margin_vec=as.vector(price_vec-cost_vec)
  
  FOC= as.matrix(ownership_strip*logit_jacobian.mean)%*%margin_vec+logit_probas.mean
  #Results:
  logit_FOC_NE= FOC
  return(logit_FOC_NE)
}

logit_FOC_NE_markup= function(beta_draws,X,price_vec,cost_vec,position,ownership,prod_names) {
  
  #Indicate inside goods
  inside_row=which(rowSums(ownership) != 0, arr.ind = TRUE)
  #update prices, only inside goods
  X$price[inside_row]=price_vec[inside_row]
  
  
  logit_probas.mean = logit_probas(beta_draws,X,prod_names)$mean
  logit_jacobian.mean = logit_jacobian(beta_draws,X,position,prod_names)$mean
  
  
  
  #strip outside good from matrices
  strip_col=which(colSums(ownership) == 0, arr.ind = TRUE)
  strip_row=which(rowSums(ownership) == 0, arr.ind = TRUE)
  logit_jacobian.mean= logit_jacobian.mean[-strip_row,-strip_col]
  ownership_strip= ownership[-strip_row,-strip_col]
  logit_probas.mean=  logit_probas.mean[-strip_row]
  price_vec=price_vec[-strip_row]
  MC=cost_vec
  
  FOC=price_vec-MC-(-solve(ownership_strip*logit_jacobian.mean)%*%logit_probas.mean) ### Similar as Eq. (6) in Morrow and Skerlos, 2011
  #Results:
  logit_FOC_NE_markup=FOC
  return(logit_FOC_NE_markup)
  
}



MC_logit_NE= function(beta_draws,X,price_vec,position,ownership,prod_names,firm_names) {
  
  #update price
  X$price=price_vec
  
  logit_probas.mean = logit_probas(beta_draws,X,prod_names)$mean
  logit_jacobian.mean = logit_jacobian(beta_draws,X,position,prod_names)$mean
  
  
  
  #strip outside good from matrices
  strip_col=which(colSums(ownership) == 0, arr.ind = TRUE)
  strip_row=which(rowSums(ownership) == 0, arr.ind = TRUE)
  logit_jacobian.mean= logit_jacobian.mean[-strip_row,-strip_col]
  ownership_strip= ownership[-strip_row,-strip_col]
  logit_probas.mean=  logit_probas.mean[-strip_row]
  price_vec=price_vec[-strip_row]
  
  
  MC=price_vec-(-solve(ownership_strip*logit_jacobian.mean)%*%logit_probas.mean)
  Margin=price_vec-MC
  #Results:
  MC_logit_NE <- NULL
  MC_logit_NE$MC = MC
  MC_logit_NE$Margin = Margin
  MC_logit_NE$products=prod_names[-strip_row]
  MC_logit_NE$firms=firm_names[-strip_row]
  return(MC_logit_NE)
  
}



