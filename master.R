setwd("~/git/BGsparse")
set.seed(123)

library(glasso)
library(statmod)
library(Kmisc)
library(Rcpp)
library(MASS)
library(caTools)
sourceCpp("colMin_Max.cpp", rebuild=TRUE, verbose=FALSE)
sourceCpp("sample_gamma.cpp", rebuild=TRUE, verbose=FALSE)
sourceCpp("sample_gamma_exact.cpp", rebuild=TRUE, verbose=FALSE)
sourceCpp("sample_Omega.cpp", rebuild=TRUE, verbose=FALSE)


source("init.R")

system.time(for (tt in 1:(burnin+niter)) {  
  if(tt%%1000==0){
    cat("iteration = ", iter <- tt, "\n")   
  }

  #### sample gamma (together with q) ####
  res1 <- sample_gamma(Omega_star_inv,Omega_star, S/n,lambda,n,gamma,Omega_star[upperind],p_gamma,pg1,pg2,r_bar,q,a_q,b_q,indi)
  gamma <- res1$gamma_tt  
  q <- res1$q
  a_gm <- res1$accept
  
  ## Gibbs sampling for gamma when p is small (enumerate all possible combination of gamma vector)
#   res1 <- sample_gamma_exact(Omega_star_inv, Omega_star, S/n, lambda, n, gamma, Omega_star[upperind],q, a_q, b_q, indi, d_sub) 
#   gamma <- res1$gamma_tt ##gamma <-res1$gamma_max (get MAP estimate)
#   q <- res1$q
  
  gamma_M <-diag(rep(1,p)); gamma_M[upperind] <- gamma; gamma_M <- gamma_M + t(gamma_M) - diag(diag(gamma_M))
  
  #### sample Omega (together with lambda) ####
  res2 <- sample_Omega(S, ind_noi_all,lambda, n, gamma_M, tau, Sigma, Omega, apost, b_lambda)
  Omega <- res2$Omega
  Sigma <- res2$Sigma
 # lambda <- res2$lambda
  
  #### sample tau ####
  lambda_tau <- lambda^2
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      if(gamma_M[i,j] == 1){
        mu_prime <- min(lambda/max(abs(Omega[i,j]), 10^-6), 10^12)
        tau[i,j] <- 1/rinvgauss(1, mu_prime, lambda_tau)
        tau[j,i] <- tau[i,j]
      }
    }
  }
  
  #### update Omega_star ###
  gla <- glasso(S/n, rho = lambda/n)
  Omega_star <- gla$wi
  Omega_star_inv <- gla$w
  if(p <= 4){
    temp <- which(Omega_star[upperind]==0)
    if( length(temp) == 0){
      d_sub <- d
    }else{
      if(length(temp) ==1){
        d_sub <- d[-which(d[,temp]>0),]
      }else{
        d_sub <- d[-which(rowSums(d[,temp])>0),]
      }
    }  
  }
  
  ######################################
  if(tt > burnin) {
    index <- tt-burnin
    gamma_sample[index,] <- gamma
    Sigma_sample[[index]] <- Sigma 
    Omega_sample[[index]] <- Omega
    lambda_sample[index] <- lambda
    #A_gamma[index]<- a_gm
  }
}
)

######### Posterior summary ############
library(plyr)
library(data.table)
Sgamma <- as.data.table(gamma_sample)
Freq <- count(Sgamma, colnames(Sgamma)) # Frequency table of sampled gamma's

row.is.a.match <- apply(gamma_sample, 1, identical, unname(unlist(Freq[which.max(Freq$freq),1:choose(p,2)])))
match.idx <- which(row.is.a.match)
total.matches <- sum(row.is.a.match)

Omega_mean <-Reduce(`+`, Omega_sample[match.idx])/total.matches
Sigma_mean <- Reduce(`+`, Sigma_sample[match.idx])/total.matches
Sigma_mean
