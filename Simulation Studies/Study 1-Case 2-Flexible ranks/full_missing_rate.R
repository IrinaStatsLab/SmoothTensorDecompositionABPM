library(Rcpp)
library(RcppArmadillo)
library(doParallel)
library(doRNG)
library(rTensor)
library(MASS)
library(refund)
library(SmoothHOOI)

source("./fpca_run.R") # script for running univariate FPCA
source("./cp_run.R")   # script for running CP decomposition with imputation and smoothness
source("./mfpca_run.R") # script for running MFPCA

## Should run this on HPC 

nworkers <- detectCores()
cl <- makePSOCKcluster(nworkers)
registerDoParallel(cl)


load("./synthetic_raw.Rda") 

N <- 100
p <- 200
D2 <- SecDiffMat(24)
rank_grid <- as.matrix(expand.grid(r1<-seq(2,6,by=1), r2<-c(2,3)))
lambda_seq <- seq(0,50,by=1)


## no missingness
set.seed(123456)

full_miss0 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","multiway","funData","MFPCA","SmoothHOOI")) %dorng%{
    sim_data <- sim_data2(L = L_tilde, R= R_tilde, mean_G = mean_G, cov_G = cov_G, 
                          E = E, p=200, noise_level=1, pattern="random", percent=0)
    Mmiss <- sim_data$sim_Mmiss
    Msmooth <- sim_data$sim_Msmooth
    
    ## Oracle
    oracle_opt <- oracle(tnsr=Mmiss@data, smooth_tnsr=Msmooth@data, rank_grid=rank_grid, lambda_seq=lambda_seq,
                         init=0, D = D2, max_iter = 500, tol = 0.1, L0 = NULL)
    
    oracle_res <- mglram(Mmiss@data, ranks=as.numeric(oracle_opt$opt_para[1:2]), lambda=as.numeric(oracle_opt$opt_para[3]),
                         L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
    
    oracle_tilde <- MakeIdent(L=oracle_res$L, G=oracle_res$G, R=oracle_res$R)
    
    oracle_loss <- loss(tnsr = oracle_res$est, smooth_tnsr=Msmooth@data)
    
    ## 10-fold cross-validation
    kcv_opt <- kcv(tnsr=Mmiss@data, rank_grid=rank_grid, lambda_seq=lambda_seq, k=10,
                   L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
    
    kcv_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_opt$opt_para[1:2]), lambda=as.numeric(kcv_opt$opt_para[3]),
                      L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
    
    kcv_tilde <- MakeIdent(L=kcv_res$L, G=kcv_res$G, R=kcv_res$R)
    
    kcv_loss <- loss(tnsr = kcv_res$est, smooth_tnsr=Msmooth@data)
    
    ## FPCA, allowing the number of principal components to vary
    fpca_res <- fpca_run(Mmiss)
    
    fpca_loss <- loss(tnsr = fpca_res$est, smooth_tnsr=Msmooth@data)
    
    ## CP with imputation and smoothness and orthogonality
    cp_res <- cp_run(Mmiss@data, R_seq=seq(2,6,by=1), const=c("ortsmo", "uncons", "uncons"))

    cp_loss <- loss(tnsr = cp_res$est, smooth_tnsr=Msmooth@data)

    ## CP with imputation and smoothness, no orthogonality
    cp_no_res <- cp_run(Mmiss@data, R_seq=seq(2,6,by=1), const=c("smooth", "uncons", "uncons"))

    cp_no_loss <- loss(tnsr = cp_no_res$est, smooth_tnsr=Msmooth@data)

    ## MFPCA
    mfpca_res <- mfpca_run(Mmiss, M_seq=seq(2,6,by=1))

    mfpca_loss <- loss(tnsr = mfpca_res$est, smooth_tnsr=Msmooth@data)
    
    
    output <- list("oracle_rank"=oracle_opt$opt_para[1], "oracle_loss"=oracle_loss,
                   "kcv_rank"=kcv_opt$opt_para[1], "kcv_loss"=kcv_loss, 
                   "fpca_rank"=ncol(fpca_res$Lfpca), "fpca_loss"=fpca_loss,
                   "cp_rank"=cp_res$best_R, "cp_loss"=cp_loss,
                   "cp_no_rank"=cp_no_res$best_R, "cp_no_loss"=cp_no_loss,
                   "mfpca_rank"=ncol(mfpca_res$Lmfpca), "mfpca_loss"=mfpca_loss
                   )
    output
  
}

save(full_miss0, file="./full_new_v2/full_miss0.Rda")


## missing rate 10%
set.seed(234567)

full_miss10 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","multiway","funData","MFPCA","SmoothHOOI")) %dorng%{
    sim_data <- sim_data2(L = L_tilde, R= R_tilde, mean_G = mean_G, cov_G = cov_G,
                          E = E, p=200, noise_level=1, pattern="random", percent=0.1)
    Mmiss <- sim_data$sim_Mmiss
    Msmooth <- sim_data$sim_Msmooth
  
    oracle_opt <- oracle(tnsr=Mmiss@data, smooth_tnsr=Msmooth@data, rank_grid=rank_grid, lambda_seq=lambda_seq,
                         init=0, D = D2, max_iter = 500, tol = 0.1, L0 = NULL)
  
    oracle_res <- mglram(Mmiss@data, ranks=as.numeric(oracle_opt$opt_para[1:2]), lambda=as.numeric(oracle_opt$opt_para[3]),
                         L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
    oracle_tilde <- MakeIdent(L=oracle_res$L, G=oracle_res$G, R=oracle_res$R)
  
    oracle_loss <- loss(tnsr = oracle_res$est, smooth_tnsr=Msmooth@data)
  
    kcv_opt <- kcv(tnsr=Mmiss@data, rank_grid=rank_grid, lambda_seq=lambda_seq, k=10,
                   L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
    kcv_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_opt$opt_para[1:2]), lambda=as.numeric(kcv_opt$opt_para[3]),
                      L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
    kcv_tilde <- MakeIdent(L=kcv_res$L, G=kcv_res$G, R=kcv_res$R)
  
    kcv_loss <- loss(tnsr = kcv_res$est, smooth_tnsr=Msmooth@data)
  
    ## FPCA, allowing the number of principal components to vary
    fpca_res <- fpca_run(Mmiss)
    
    fpca_loss <- loss(tnsr = fpca_res$est, smooth_tnsr=Msmooth@data)
    
    ## CP with imputation and smoothness and orthogonality
    cp_res <- cp_run(Mmiss@data, R_seq=seq(2,6,by=1), const=c("ortsmo", "uncons", "uncons"))

    cp_loss <- loss(tnsr = cp_res$est, smooth_tnsr=Msmooth@data)

    ## CP with imputation and smoothness, no orthogonality
    cp_no_res <- cp_run(Mmiss@data, R_seq=seq(2,6,by=1), const=c("smooth", "uncons", "uncons"))

    cp_no_loss <- loss(tnsr = cp_no_res$est, smooth_tnsr=Msmooth@data)

    ## MFPCA
    mfpca_res <- mfpca_run(Mmiss, M_seq=seq(2,6,by=1))

    mfpca_loss <- loss(tnsr = mfpca_res$est, smooth_tnsr=Msmooth@data)
    
    
    output <- list("oracle_rank"=oracle_opt$opt_para[1], "oracle_loss"=oracle_loss,
                   "kcv_rank"=kcv_opt$opt_para[1], "kcv_loss"=kcv_loss, 
                   "fpca_rank"=ncol(fpca_res$Lfpca), "fpca_loss"=fpca_loss,
                   "cp_rank"=cp_res$best_R, "cp_loss"=cp_loss,
                   "cp_no_rank"=cp_no_res$best_R, "cp_no_loss"=cp_no_loss,
                   "mfpca_rank"=ncol(mfpca_res$Lmfpca), "mfpca_loss"=mfpca_loss
                   )
    output
}

save(full_miss10, file="./full_new_v2/full_miss10.Rda")


## missing rate 20%
set.seed(345678)

full_miss20 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","multiway","funData","MFPCA","SmoothHOOI")) %dorng%{
    sim_data <- sim_data2(L = L_tilde, R= R_tilde, mean_G = mean_G, cov_G = cov_G,
                          E = E, p=200, noise_level=1, pattern="random", percent=0.2)
    Mmiss <- sim_data$sim_Mmiss
    Msmooth <- sim_data$sim_Msmooth

    oracle_opt <- oracle(tnsr=Mmiss@data, smooth_tnsr=Msmooth@data, rank_grid=rank_grid, lambda_seq=lambda_seq,
                         init=0, D = D2, max_iter = 500, tol = 0.1, L0 = NULL)

    oracle_res <- mglram(Mmiss@data, ranks=as.numeric(oracle_opt$opt_para[1:2]), lambda=as.numeric(oracle_opt$opt_para[3]),
                         L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)

    oracle_tilde <- MakeIdent(L=oracle_res$L, G=oracle_res$G, R=oracle_res$R)

    oracle_loss <- loss(tnsr = oracle_res$est, smooth_tnsr=Msmooth@data)

    kcv_opt <- kcv(tnsr=Mmiss@data, rank_grid=rank_grid, lambda_seq=lambda_seq, k=10,
                   L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)

    kcv_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_opt$opt_para[1:2]), lambda=as.numeric(kcv_opt$opt_para[3]),
                      L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)

    kcv_tilde <- MakeIdent(L=kcv_res$L, G=kcv_res$G, R=kcv_res$R)

    kcv_loss <- loss(tnsr = kcv_res$est, smooth_tnsr=Msmooth@data)

    ## FPCA, allowing the number of principal components to vary
    fpca_res <- fpca_run(Mmiss)

    fpca_loss <- loss(tnsr = fpca_res$est, smooth_tnsr=Msmooth@data)

    ## CP with imputation and smoothness and orthogonality
    cp_res <- cp_run(Mmiss@data, R_seq=seq(2,6,by=1), const=c("ortsmo", "uncons", "uncons"))

    cp_loss <- loss(tnsr = cp_res$est, smooth_tnsr=Msmooth@data)

    ## CP with imputation and smoothness, no orthogonality
    cp_no_res <- cp_run(Mmiss@data, R_seq=seq(2,6,by=1), const=c("smooth", "uncons", "uncons"))

    cp_no_loss <- loss(tnsr = cp_no_res$est, smooth_tnsr=Msmooth@data)

    ## MFPCA
    mfpca_res <- mfpca_run(Mmiss, M_seq=seq(2,6,by=1))

    mfpca_loss <- loss(tnsr = mfpca_res$est, smooth_tnsr=Msmooth@data)

    output <- list("oracle_rank"=oracle_opt$opt_para[1], "oracle_loss"=oracle_loss,
                   "kcv_rank"=kcv_opt$opt_para[1], "kcv_loss"=kcv_loss,
                   "fpca_rank"=ncol(fpca_res$Lfpca), "fpca_loss"=fpca_loss,
                   "cp_rank"=cp_res$best_R, "cp_loss"=cp_loss,
                   "cp_no_rank"=cp_no_res$best_R, "cp_no_loss"=cp_no_loss,
                   "mfpca_rank"=ncol(mfpca_res$Lmfpca), "mfpca_loss"=mfpca_loss
    )
    output
}

save(full_miss20, file="./full_new_v2/full_miss20.Rda")


## missing rate 50%
set.seed(456789)

full_miss50 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","multiway","funData","MFPCA","SmoothHOOI")) %dorng%{
    sim_data <- sim_data2(L = L_tilde, R= R_tilde, mean_G = mean_G, cov_G = cov_G,
                          E = E, p=200, noise_level=1, pattern="random", percent=0.5)
    Mmiss <- sim_data$sim_Mmiss
    Msmooth <- sim_data$sim_Msmooth

    oracle_opt <- oracle(tnsr=Mmiss@data, smooth_tnsr=Msmooth@data, rank_grid=rank_grid, lambda_seq=lambda_seq,
                         init=0, D = D2, max_iter = 500, tol = 0.1, L0 = NULL)

    oracle_res <- mglram(Mmiss@data, ranks=as.numeric(oracle_opt$opt_para[1:2]), lambda=as.numeric(oracle_opt$opt_para[3]),
                         L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)

    oracle_tilde <- MakeIdent(L=oracle_res$L, G=oracle_res$G, R=oracle_res$R)

    oracle_loss <- loss(tnsr = oracle_res$est, smooth_tnsr=Msmooth@data)

    kcv_opt <- kcv(tnsr=Mmiss@data, rank_grid=rank_grid, lambda_seq=lambda_seq, k=10,
                   L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)

    kcv_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_opt$opt_para[1:2]), lambda=as.numeric(kcv_opt$opt_para[3]),
                      L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)

    kcv_tilde <- MakeIdent(L=kcv_res$L, G=kcv_res$G, R=kcv_res$R)

    kcv_loss <- loss(tnsr = kcv_res$est, smooth_tnsr=Msmooth@data)

    ## FPCA, allowing the number of principal components to vary
    fpca_res <- fpca_run(Mmiss)

    fpca_loss <- loss(tnsr = fpca_res$est, smooth_tnsr=Msmooth@data)

    ## CP with imputation and smoothness and orthogonality
    cp_res <- cp_run(Mmiss@data, R_seq=seq(2,6,by=1), const=c("ortsmo", "uncons", "uncons"))

    cp_loss <- loss(tnsr = cp_res$est, smooth_tnsr=Msmooth@data)

    ## CP with imputation and smoothness, no orthogonality
    cp_no_res <- cp_run(Mmiss@data, R_seq=seq(2,6,by=1), const=c("smooth", "uncons", "uncons"))

    cp_no_loss <- loss(tnsr = cp_no_res$est, smooth_tnsr=Msmooth@data)

    ## MFPCA
    mfpca_res <- mfpca_run(Mmiss, M_seq=seq(2,6,by=1))

    mfpca_loss <- loss(tnsr = mfpca_res$est, smooth_tnsr=Msmooth@data)

    output <- list("oracle_rank"=oracle_opt$opt_para[1], "oracle_loss"=oracle_loss,
                   "kcv_rank"=kcv_opt$opt_para[1], "kcv_loss"=kcv_loss,
                   "fpca_rank"=ncol(fpca_res$Lfpca), "fpca_loss"=fpca_loss,
                   "cp_rank"=cp_res$best_R, "cp_loss"=cp_loss,
                   "cp_no_rank"=cp_no_res$best_R, "cp_no_loss"=cp_no_loss,
                   "mfpca_rank"=ncol(mfpca_res$Lmfpca), "mfpca_loss"=mfpca_loss
    )
    output
}

save(full_miss50, file="./full_new_v2/full_miss50.Rda")

stopCluster(cl)
