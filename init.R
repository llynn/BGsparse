p <- 30; n <- 100; # p: dimension, n: sample size
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
CTrue[abs(CTrue) < 1.0e-10] = 0
#### generate data
Y <- t(chol(SigTrue))%*%matrix(rnorm(p*n), nrow = p) # Y_i ~ N(0, SigTrue)
S <- Y%*%t(Y)

#### initialization parameters ####
a_lambda <- 1; b_lambda <- 0.1; # hyperparameters for lambda ~ Ga(a_lambda, b_lambda)
Sigma_sample <- vector(mode="list", length = niter)
Sigma <- S/n


Omega_sample <- Sigma_sample
Omega <- solve(Sigma)

lambda_sample <- rep(0,niter) #parameter for both Exp and DE for \omega_ij and \omega_ii 
apost <- a_lambda + p*(p+1)/2
bpost <- b_lambda + sum(abs(Omega))/2
#lambda <- rgamma(1,apost,1/bpost)
lambda = 1

#############################################
gamma_sample <- array(0, dim = c(niter, choose(p,2))) #choose(p,2) edge inclusion indicators (excluding diagonal elements)
gamma <- rep(as.integer(1), choose(p,2))
names(gamma) <- upperind


Omega_star <- glasso(S/n, rho = lambda/n)$wi # graphical lasso solution
Omega_star_inv <- chol2inv(chol(Omega_star))

###fixed_position: which position gamma == 0
fixed_position <- which(Omega_star[upperind] == 0)
freemove <- 1:choose(p,2)

if(length(fixed_position)>0){
  freemove <- freemove[-fixed_position]
  gamma[fixed_position] = 0
}

gamma_M <-diag(rep(1,p))
gamma_M[upperind] <- gamma

gamma_M <- gamma_M + t(gamma_M) - diag(diag(gamma_M)) 
p_gamma <- 0.5 # parameters for propose new gamma
pg1 = 0.2;
pg2 = 0.6;

q <- 0.4
a_q <- .4; b_q <- .5; # parameter for distribution of q : q ~ beta(a_q, b_q)1_{q < 0.5}

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





