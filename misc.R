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


/home/llin/tophat-2.0.11.Linux_x86_64/tophat2 -p 1 -G lncRNA.e67xloc.new.gtf --keep-fasta-order --no-novel-juncs /shared/silo_researcher/Gottardo_R/jingyuan_working/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/BowtieIndex/genome SRR073723_1.fastq

