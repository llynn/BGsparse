#### sample gamma only p<=4 #########
source("init.R")

res1 <- sample_gamma_exact(Omega_star_inv, Omega_star, S/n, 0.1*n, n, Omega_star[upperind],0.4, a_q, b_q, indi, d_sub) 

#d_sub: all possible structures, HD: hamming distance to glasso
HD <- apply(d_sub, 1, function (x) {
  sum(x != gamma_lasso)
})

median(res1$weights[which(HD <=2)]) # median weights for models within a HD of 2 from glasso
median(res1$weights)

which.d_sub <- apply(d_sub, 1, identical, Tgamma)
res1$weights[which(which.d_sub)]  #weights for the true structure (gamma)


