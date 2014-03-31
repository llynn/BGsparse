setwd("~/git/BGsparse")
set.seed(123)

library(glasso)
library(statmod)
library(Kmisc)
library(Rcpp)


sourceCpp("colMin_Max.cpp", rebuild=TRUE, verbose=FALSE)
sourceCpp("sample_gamma.cpp", rebuild=TRUE, verbose=FALSE)
sourceCpp("sample_Omega.cpp", rebuild=TRUE, verbose=FALSE)

burnin  = 1000; niter= 2000; 
source("init.R")

for (tt in 1:(burnin+niter)) {  
  #cat("iteration = ", iter <- tt, "\n")
  #### sample gamma (together with q) ####

  system.time(res1 <- sample_gamma(Omega_star_inv,S,lambda,n,gamma,Omega_star[upperind],p_gamma,pg1,pg2,r_bar,q,a_q,b_q,indi))
  gamma <- res1$gamma_tt
  q <- res1$q
  a_gm <- res1$accept
  gamma_M <-diag(rep(1,p)); gamma_M[upperind] <- gamma; gamma_M <- gamma_M + t(gamma_M) - diag(diag(gamma_M))
  
  
  #### sample Omega (together with lambda) ####
  lambda<- sample_Omega(Omega_star, S, ind_noi_all,lambda,n, gamma_M, Omega_star[upperind],tau,Sigma, indi, Omega, apost, b_lambda)

  #### sample tau ####
  lambda_tau <- lambda^2
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      if(gamma_M[i,j] == 1){
        mu_prime <- min(lambda/max(abs(Omega[i,j])), 10^-6)
        tau[i,j] <- 1/rinvgauss(1, mu_prime, lambda_tau)
        tau[j,i] <- tau[i,j]
      }
#       else{
#         tau[i,j] <- rexp(1, lambda_tau/2)
#         tau[j,i] <- tau[i,j]
#       }
    }
  }
  
  #### update Omega_star 
  Omega_star <- glasso(S/n, rho = lambda/n)$wi
  Omega_star_inv <- chol2inv(chol(Omega_star))
  ######################################
  if(tt > burnin) {
    index <- tt-burnin
    gamma_sample[index,] <- gamma
    Sigma_sample[[index]] <- Sigma 
    Omega_sample[[index]] <- Omega
    lambda_sample[index] <- lambda
    A_gamma[index]<- a_gm
  }
}
