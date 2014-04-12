library(rbenchmark)

n <- 1e5


system.time(res <- benchmark(determinant(Omega_star),
                 R_det_(Omega_star),
                 replications = 50000))
print(res[,1:4])


system.time(output<-det_H_omega(Omega_star_inv, gamma, sum(gamma)+p, p, indi, choose(p,2)))


system.time(output_slow <-det_H_omega_slow(Omega_star_inv, gamma, sum(gamma)+p, p, indi, choose(p,2)))


res<- sample_Omega(S, ind_noi_all,lambda, n, gamma_M, tau, Sigma, Omega, apost, b_lambda)

res_W<- sample_Omega_W(S, ind_noi_all,lambda, n, gamma_M, tau, Sigma, Omega, apost, b_lambda)

