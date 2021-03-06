burnin  = 10000; niter= 10000;  
p <- 4; n <- 100; # p: dimension, n: sample size
r_bar <- choose(p,2) #loose bound first, fixed

indmx <- matrix(1:p^2,p,p) 
upperind <- indmx[upper.tri(indmx)] 
lowerind <- indmx[lower.tri(indmx)] 

indi <- data.frame(which(upper.tri(indmx, diag = T), arr.ind = T))
indi <- as.matrix(indi[with(indi, order(row)),]-1)

ind_noi_all <- matrix(0, nrow = p-1, ncol = p)
for (i in 1:p){
  if (i == 1){
    ind_noi <- t(2:p)
  }else if (i == p){
    ind_noi <- t(1:(p-1))
  }else{
    ind_noi = t(c(1:(i-1),(i+1):p))
  }
  ind_noi_all[,i] <- ind_noi
}
ind_noi_all <- ind_noi_all-1

#### AR(1) case
SigTrue <- toeplitz(0.7^(0:(p-1)))
CTrue <- chol2inv(chol(SigTrue))
CTrue[abs(CTrue) <= 10^-10] = 0
Tgamma <- matrix(1L, nrow(CTrue), ncol(CTrue))
Tgamma[CTrue ==0] = 0L
Tgamma = Tgamma[upperind]
#### AR(2) case
# CTrue  <- toeplitz(c(1,0.5,0.25,rep(0,p-3)))
# SigTrue = chol2inv(chol(CTrue))

#### generate data
#Y <- t(chol(SigTrue))%*%matrix(rnorm(p*n), nrow = p) # Y_i ~ N(0, SigTrue)
Y <- t(mvrnorm(n, rep(0,p), SigTrue))
S <- Y%*%t(Y)

#### initialization parameters ####
a_lambda <- 1; b_lambda <- .1; # hyperparameters for lambda ~ Ga(a_lambda, b_lambda)
Sigma_sample <- vector(mode="list", length = niter)
Sigma <- S/n


Omega_sample <- Sigma_sample
Omega <- solve(Sigma)

#############################################
gamma_sample <- array(0, dim = c(niter, choose(p,2))) #choose(p,2) edge inclusion indicators (excluding diagonal elements)
gamma <- rep(as.integer(1), choose(p,2))
names(gamma) <- upperind

lambda_sample <- rep(0,niter) #parameter for both Exp and DE for \omega_ij and \omega_ii 
#apost <- a_lambda + p*(p+1)/2
apost <- a_lambda + sum(gamma)+p
bpost <- b_lambda + sum(abs(Omega))/2
#lambda <- rgamma(1,apost,1/bpost)
lambda = 1


#Omega_star <- glasso(S/n, rho = lambda/n)$wi # graphical lasso solution for precision matrix
Omega_star <- glasso(S/n, rho = .1)$wi
Omega_star_inv <- chol2inv(chol(Omega_star))

###fixed_position: which position gamma == 0
fixed_position <- which(Omega_star[upperind] == 0)
freemove <- 1:choose(p,2)

if(length(fixed_position)>0){
  freemove <- freemove[-fixed_position]
  gamma[fixed_position] = 0
}

gamma_lasso <- gamma

gamma_M <-diag(rep(1,p))
gamma_M[upperind] <- gamma

gamma_M <- gamma_M + t(gamma_M) - diag(diag(gamma_M)) 

p_gamma <- 0.5 # parameters for propose new gamma
pg1 = 0.1;
pg2 = 0.4;

q <- 0.5
#a_q <- (length(freemove) - length(fixed_position)); b_q <- length(freemove); # parameter for distribution of q : q ~ beta(a_q, b_q)1_{q < 0.5}
a_q <- 2; b_q <- 5;
##### sample tau off-diagonal based on gamma
tau <- array(0, dim = c(p,p))
Oadjust <- colMin_Max(rbind(abs(Omega[upperind]),10^-6),1)   
mu_tau <- colMin_Max(rbind(lambda/Oadjust,10^12), 0)
lambda_tau <- lambda^2
tau_temp <- 1/unlist(lapply(mu_tau, function(x) rinvgauss(1, x, lambda_tau)))
tau[upperind] <- tau_temp
tau[lowerind] <- tau_temp

if(length(fixed_position) > 0){
  tmp <- rexp(length(fixed_position), lambda^2/2)
  tau[upperind[fixed_position]] = tmp
  tau[lowerind[fixed_position]] = tmp
}

####### acceptance rate ###########
A_gamma <- rep(0, niter)

###########################
if(p <=4) {
  K <- 2^choose(p,2)
  M <- choose(p,2)
  d = array(as.integer(0),dim=c(K, M))
  count=1;
  for (kk in 1:M) {
    set = combs(1:M,kk);
    ll = dim(combs(1:M,kk)); 
    
    for ( lll in 1:ll[1]) {
      d[(count-1+lll),set[lll,]] = as.integer(1);
    }
    
    count = count+choose(M,kk);
  }
  temp <- which(rowSums(d) <= r_bar)
  d <- d[temp,]
  
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





