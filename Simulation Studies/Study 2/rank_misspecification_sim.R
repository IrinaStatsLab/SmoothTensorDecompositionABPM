library(Rcpp)
library(RcppArmadillo)
library(doParallel)
library(doRNG)
library(rTensor)
library(MASS)
library(refund)
library(SmoothHOOI)
library(MFPCA)
library(multiway)

source("./rank_truncation.R")
source("./calc_EV.R")


### Should run this on HPC 

nworkers <- detectCores()
cl <- makePSOCKcluster(nworkers)
registerDoParallel(cl)

# Results generated from running SmoothHOOI algorithm on real data, with hyperparameter optimization and identifiability correction
load("./synthetic_raw.Rda") 

N <- 100 # Number of replicates 
D2 <- SecDiffMat(24) # Second-order difference matrix

# true rank is (3,2)
true_rank_grid <- as.matrix(expand.grid(r1<-c(3), r2<-c(2)))
# other rank combinations
rank_grid_42 <- as.matrix(expand.grid(r1<-c(4), r2<-c(2)))
rank_grid_62 <- as.matrix(expand.grid(r1<-c(6), r2<-c(2)))
rank_grid_33 <- as.matrix(expand.grid(r1<-c(3), r2<-c(3)))
rank_grid_43 <- as.matrix(expand.grid(r1<-c(4), r2<-c(3)))
rank_grid_63 <- as.matrix(expand.grid(r1<-c(6), r2<-c(3)))

# lambda range
lambda_seq <- seq(0,50,by=1)

###################################################################

# no missingness
set.seed(20261001)

miss0 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","multiway","funData","MFPCA","SmoothHOOI")) %dorng%{
  ## Generate synthetic data
  sim_data <- sim_data2(L = L_tilde, R= R_tilde, mean_G = mean_G, cov_G = cov_G,
                        E = E, p=200, noise_level=1, pattern="random", percent=0)
  Mmiss <- sim_data$sim_Mmiss
  Msmooth <- sim_data$sim_Msmooth
  
  ## 10-fold cross-validation--true rank
  kcv_true_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=true_rank_grid, lambda_seq=lambda_seq, k=10,
                             L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_true_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_true_opt$opt_para[1:2]), lambda=as.numeric(kcv_true_opt$opt_para[3]),
                         L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_true_tilde <- MakeIdent(L=kcv_true_res$L, G=kcv_true_res$G, R=kcv_true_res$R)
  
  kcv_true_loss <- loss(tnsr = kcv_true_res$est, true_tnsr=Msmooth@data,
                        L = kcv_true_tilde$L_tilde, true_L = L_tilde,
                        R = kcv_true_tilde$R_tilde, true_R = R_tilde)
  
  EV_true <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_true_tilde$L_tilde,
    G_tilde = kcv_true_tilde$G_tilde,
    R_tilde = kcv_true_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--42
  kcv_42_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_42, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_42_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_42_opt$opt_para[1:2]), lambda=as.numeric(kcv_42_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_42_tilde <- MakeIdent(L=kcv_42_res$L, G=kcv_42_res$G, R=kcv_42_res$R)
  
  kcv_42_lossM <- loss(tnsr = kcv_42_res$est, true_tnsr=Msmooth@data)
  
  kcv_42_truncated_res <- rank_truncation(Mmiss, mode="L", true_rank_L=3, true_rank_R=NULL,
                                          L_tilde=kcv_42_tilde$L_tilde, G_tilde=kcv_42_tilde$G_tilde, R_tilde=kcv_42_tilde$R_tilde)
  
  kcv_42_loss_truncated <- loss(tnsr = kcv_42_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_42_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_42_tilde$R_tilde, true_R = R_tilde)
  
  EV_42 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_42_tilde$L_tilde,
    G_tilde = kcv_42_tilde$G_tilde,
    R_tilde = kcv_42_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--62
  kcv_62_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_62, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_62_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_62_opt$opt_para[1:2]), lambda=as.numeric(kcv_62_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_62_tilde <- MakeIdent(L=kcv_62_res$L, G=kcv_62_res$G, R=kcv_62_res$R)
  
  kcv_62_lossM <- loss(tnsr = kcv_62_res$est, true_tnsr=Msmooth@data)
  
  kcv_62_truncated_res <- rank_truncation(Mmiss, mode="L", true_rank_L=3, true_rank_R=NULL,
                                          L_tilde=kcv_62_tilde$L_tilde, G_tilde=kcv_62_tilde$G_tilde, R_tilde=kcv_62_tilde$R_tilde)
  
  kcv_62_loss_truncated <- loss(tnsr = kcv_62_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_62_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_62_tilde$R_tilde, true_R = R_tilde)
  
  EV_62 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_62_tilde$L_tilde,
    G_tilde = kcv_62_tilde$G_tilde,
    R_tilde = kcv_62_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--33
  kcv_33_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_33, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_33_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_33_opt$opt_para[1:2]), lambda=as.numeric(kcv_33_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_33_tilde <- MakeIdent(L=kcv_33_res$L, G=kcv_33_res$G, R=kcv_33_res$R)
  
  kcv_33_lossM <- loss(tnsr = kcv_33_res$est, true_tnsr=Msmooth@data)
  
  kcv_33_truncated_res <- rank_truncation(Mmiss, mode="R", true_rank_L=NULL, true_rank_R=2,
                                          L_tilde=kcv_33_tilde$L_tilde, G_tilde=kcv_33_tilde$G_tilde, R_tilde=kcv_33_tilde$R_tilde)
  
  kcv_33_loss_truncated <- loss(tnsr = kcv_33_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_33_tilde$L_tilde, true_L = L_tilde,
                                R = kcv_33_truncated_res$R_truncated, true_R = R_tilde)
  
  EV_33 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_33_tilde$L_tilde,
    G_tilde = kcv_33_tilde$G_tilde,
    R_tilde = kcv_33_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--43
  kcv_43_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_43, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_43_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_43_opt$opt_para[1:2]), lambda=as.numeric(kcv_43_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_43_tilde <- MakeIdent(L=kcv_43_res$L, G=kcv_43_res$G, R=kcv_43_res$R)
  
  kcv_43_lossM <- loss(tnsr = kcv_43_res$est, true_tnsr=Msmooth@data)
  
  kcv_43_truncated_res <- rank_truncation(Mmiss, mode="both", true_rank_L=3, true_rank_R=2,
                                          L_tilde=kcv_43_tilde$L_tilde, G_tilde=kcv_43_tilde$G_tilde, R_tilde=kcv_43_tilde$R_tilde)
  
  kcv_43_loss_truncated <- loss(tnsr = kcv_43_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_43_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_43_truncated_res$R_truncated, true_R = R_tilde)
  
  EV_43 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_43_tilde$L_tilde,
    G_tilde = kcv_43_tilde$G_tilde,
    R_tilde = kcv_43_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--63
  kcv_63_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_63, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_63_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_63_opt$opt_para[1:2]), lambda=as.numeric(kcv_63_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_63_tilde <- MakeIdent(L=kcv_63_res$L, G=kcv_63_res$G, R=kcv_63_res$R)
  
  kcv_63_lossM <- loss(tnsr = kcv_63_res$est, true_tnsr=Msmooth@data)
  
  kcv_63_truncated_res <- rank_truncation(Mmiss, mode="both", true_rank_L=3, true_rank_R=2,
                                          L_tilde=kcv_63_tilde$L_tilde, G_tilde=kcv_63_tilde$G_tilde, R_tilde=kcv_63_tilde$R_tilde)
  
  kcv_63_loss_truncated <- loss(tnsr = kcv_63_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_63_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_63_truncated_res$R_truncated, true_R = R_tilde)
  
  
  EV_63 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_63_tilde$L_tilde,
    G_tilde = kcv_63_tilde$G_tilde,
    R_tilde = kcv_63_tilde$R_tilde
  )
  
  output <- list("kcv_true_loss"=kcv_true_loss, lambda_true = as.numeric(kcv_true_opt$opt_para[3]),
                 "kcv_42_lossM"=kcv_42_lossM, "kcv_42_loss_truncated"=kcv_42_loss_truncated,"lambda_42" = as.numeric(kcv_42_opt$opt_para[3]),
                 "kcv_62_lossM"=kcv_62_lossM, "kcv_62_loss_truncated"=kcv_62_loss_truncated,"lambda_62" = as.numeric(kcv_62_opt$opt_para[3]),
                 "kcv_33_lossM"=kcv_33_lossM, "kcv_33_loss_truncated"=kcv_33_loss_truncated,"lambda_33" = as.numeric(kcv_33_opt$opt_para[3]),
                 "kcv_43_lossM"=kcv_43_lossM, "kcv_43_loss_truncated"=kcv_43_loss_truncated,"lambda_43" = as.numeric(kcv_43_opt$opt_para[3]),
                 "kcv_63_lossM"=kcv_63_lossM, "kcv_63_loss_truncated"=kcv_63_loss_truncated,"lambda_63" = as.numeric(kcv_63_opt$opt_para[3]),
                 "EV_true" = EV_true,"EV_42" = EV_42,"EV_62" = EV_62, "EV_33" = EV_33,"EV_43" = EV_43,"EV_63" = EV_63
  )
  output
}

save(miss0, file="./rank_misspecification_ev/miss0.Rda")

## missing rate 10%
set.seed(20261002)

miss10 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","multiway","funData","MFPCA","SmoothHOOI")) %dorng%{
  sim_data <- sim_data2(L = L_tilde, R= R_tilde, mean_G = mean_G, cov_G = cov_G,
                        E = E, p=200, noise_level=1, pattern="random", percent=0.1)
  Mmiss <- sim_data$sim_Mmiss
  Msmooth <- sim_data$sim_Msmooth
  
  ## 10-fold cross-validation--true rank
  kcv_true_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=true_rank_grid, lambda_seq=lambda_seq, k=10,
                             L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_true_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_true_opt$opt_para[1:2]), lambda=as.numeric(kcv_true_opt$opt_para[3]),
                         L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_true_tilde <- MakeIdent(L=kcv_true_res$L, G=kcv_true_res$G, R=kcv_true_res$R)
  
  kcv_true_loss <- loss(tnsr = kcv_true_res$est, true_tnsr=Msmooth@data,
                        L = kcv_true_tilde$L_tilde, true_L = L_tilde,
                        R = kcv_true_tilde$R_tilde, true_R = R_tilde)
  
  EV_true <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_true_tilde$L_tilde,
    G_tilde = kcv_true_tilde$G_tilde,
    R_tilde = kcv_true_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--42
  kcv_42_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_42, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_42_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_42_opt$opt_para[1:2]), lambda=as.numeric(kcv_42_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_42_tilde <- MakeIdent(L=kcv_42_res$L, G=kcv_42_res$G, R=kcv_42_res$R)
  
  kcv_42_lossM <- loss(tnsr = kcv_42_res$est, true_tnsr=Msmooth@data)
  
  kcv_42_truncated_res <- rank_truncation(Mmiss, mode="L", true_rank_L=3, true_rank_R=NULL,
                                          L_tilde=kcv_42_tilde$L_tilde, G_tilde=kcv_42_tilde$G_tilde, R_tilde=kcv_42_tilde$R_tilde)
  
  kcv_42_loss_truncated <- loss(tnsr = kcv_42_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_42_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_42_tilde$R_tilde, true_R = R_tilde)
  
  EV_42 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_42_tilde$L_tilde,
    G_tilde = kcv_42_tilde$G_tilde,
    R_tilde = kcv_42_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--62
  kcv_62_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_62, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_62_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_62_opt$opt_para[1:2]), lambda=as.numeric(kcv_62_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_62_tilde <- MakeIdent(L=kcv_62_res$L, G=kcv_62_res$G, R=kcv_62_res$R)
  
  kcv_62_lossM <- loss(tnsr = kcv_62_res$est, true_tnsr=Msmooth@data)
  
  kcv_62_truncated_res <- rank_truncation(Mmiss, mode="L", true_rank_L=3, true_rank_R=NULL,
                                          L_tilde=kcv_62_tilde$L_tilde, G_tilde=kcv_62_tilde$G_tilde, R_tilde=kcv_62_tilde$R_tilde)
  
  kcv_62_loss_truncated <- loss(tnsr = kcv_62_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_62_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_62_tilde$R_tilde, true_R = R_tilde)
  
  EV_62 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_62_tilde$L_tilde,
    G_tilde = kcv_62_tilde$G_tilde,
    R_tilde = kcv_62_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--33
  kcv_33_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_33, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_33_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_33_opt$opt_para[1:2]), lambda=as.numeric(kcv_33_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_33_tilde <- MakeIdent(L=kcv_33_res$L, G=kcv_33_res$G, R=kcv_33_res$R)
  
  kcv_33_lossM <- loss(tnsr = kcv_33_res$est, true_tnsr=Msmooth@data)
  
  kcv_33_truncated_res <- rank_truncation(Mmiss, mode="R", true_rank_L=NULL, true_rank_R=2,
                                          L_tilde=kcv_33_tilde$L_tilde, G_tilde=kcv_33_tilde$G_tilde, R_tilde=kcv_33_tilde$R_tilde)
  
  kcv_33_loss_truncated <- loss(tnsr = kcv_33_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_33_tilde$L_tilde, true_L = L_tilde,
                                R = kcv_33_truncated_res$R_truncated, true_R = R_tilde)
  
  EV_33 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_33_tilde$L_tilde,
    G_tilde = kcv_33_tilde$G_tilde,
    R_tilde = kcv_33_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--43
  kcv_43_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_43, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_43_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_43_opt$opt_para[1:2]), lambda=as.numeric(kcv_43_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_43_tilde <- MakeIdent(L=kcv_43_res$L, G=kcv_43_res$G, R=kcv_43_res$R)
  
  kcv_43_lossM <- loss(tnsr = kcv_43_res$est, true_tnsr=Msmooth@data)
  
  kcv_43_truncated_res <- rank_truncation(Mmiss, mode="both", true_rank_L=3, true_rank_R=2,
                                          L_tilde=kcv_43_tilde$L_tilde, G_tilde=kcv_43_tilde$G_tilde, R_tilde=kcv_43_tilde$R_tilde)
  
  kcv_43_loss_truncated <- loss(tnsr = kcv_43_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_43_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_43_truncated_res$R_truncated, true_R = R_tilde)
  
  EV_43 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_43_tilde$L_tilde,
    G_tilde = kcv_43_tilde$G_tilde,
    R_tilde = kcv_43_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--63
  kcv_63_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_63, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_63_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_63_opt$opt_para[1:2]), lambda=as.numeric(kcv_63_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_63_tilde <- MakeIdent(L=kcv_63_res$L, G=kcv_63_res$G, R=kcv_63_res$R)
  
  kcv_63_lossM <- loss(tnsr = kcv_63_res$est, true_tnsr=Msmooth@data)
  
  kcv_63_truncated_res <- rank_truncation(Mmiss, mode="both", true_rank_L=3, true_rank_R=2,
                                          L_tilde=kcv_63_tilde$L_tilde, G_tilde=kcv_63_tilde$G_tilde, R_tilde=kcv_63_tilde$R_tilde)
  
  kcv_63_loss_truncated <- loss(tnsr = kcv_63_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_63_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_63_truncated_res$R_truncated, true_R = R_tilde)
  
  
  EV_63 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_63_tilde$L_tilde,
    G_tilde = kcv_63_tilde$G_tilde,
    R_tilde = kcv_63_tilde$R_tilde
  )
  
  output <- list("kcv_true_loss"=kcv_true_loss, lambda_true = as.numeric(kcv_true_opt$opt_para[3]),
                 "kcv_42_lossM"=kcv_42_lossM, "kcv_42_loss_truncated"=kcv_42_loss_truncated,"lambda_42" = as.numeric(kcv_42_opt$opt_para[3]),
                 "kcv_62_lossM"=kcv_62_lossM, "kcv_62_loss_truncated"=kcv_62_loss_truncated,"lambda_62" = as.numeric(kcv_62_opt$opt_para[3]),
                 "kcv_33_lossM"=kcv_33_lossM, "kcv_33_loss_truncated"=kcv_33_loss_truncated,"lambda_33" = as.numeric(kcv_33_opt$opt_para[3]),
                 "kcv_43_lossM"=kcv_43_lossM, "kcv_43_loss_truncated"=kcv_43_loss_truncated,"lambda_43" = as.numeric(kcv_43_opt$opt_para[3]),
                 "kcv_63_lossM"=kcv_63_lossM, "kcv_63_loss_truncated"=kcv_63_loss_truncated,"lambda_63" = as.numeric(kcv_63_opt$opt_para[3]),
                 "EV_true" = EV_true,"EV_42" = EV_42,"EV_62" = EV_62, "EV_33" = EV_33,"EV_43" = EV_43,"EV_63" = EV_63
  )
  output
}

save(miss10, file="./rank_misspecification_ev/miss10.Rda")

## missing rate 20%
set.seed(20261003)

miss20 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","multiway","funData","MFPCA","SmoothHOOI")) %dorng%{
  sim_data <- sim_data2(L = L_tilde, R= R_tilde, mean_G = mean_G, cov_G = cov_G,
                        E = E, p=200, noise_level=1, pattern="random", percent=0.2)
  Mmiss <- sim_data$sim_Mmiss
  Msmooth <- sim_data$sim_Msmooth
  
  ## 10-fold cross-validation--true rank
  kcv_true_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=true_rank_grid, lambda_seq=lambda_seq, k=10,
                             L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_true_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_true_opt$opt_para[1:2]), lambda=as.numeric(kcv_true_opt$opt_para[3]),
                         L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_true_tilde <- MakeIdent(L=kcv_true_res$L, G=kcv_true_res$G, R=kcv_true_res$R)
  
  kcv_true_loss <- loss(tnsr = kcv_true_res$est, true_tnsr=Msmooth@data,
                        L = kcv_true_tilde$L_tilde, true_L = L_tilde,
                        R = kcv_true_tilde$R_tilde, true_R = R_tilde)
  
  
  EV_true <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_true_tilde$L_tilde,
    G_tilde = kcv_true_tilde$G_tilde,
    R_tilde = kcv_true_tilde$R_tilde
  )
  
  ## 10-fold cross-validation--42
  kcv_42_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_42, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_42_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_42_opt$opt_para[1:2]), lambda=as.numeric(kcv_42_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_42_tilde <- MakeIdent(L=kcv_42_res$L, G=kcv_42_res$G, R=kcv_42_res$R)
  
  kcv_42_lossM <- loss(tnsr = kcv_42_res$est, true_tnsr=Msmooth@data)
  
  kcv_42_truncated_res <- rank_truncation(Mmiss, mode="L", true_rank_L=3, true_rank_R=NULL,
                                          L_tilde=kcv_42_tilde$L_tilde, G_tilde=kcv_42_tilde$G_tilde, R_tilde=kcv_42_tilde$R_tilde)
  
  kcv_42_loss_truncated <- loss(tnsr = kcv_42_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_42_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_42_tilde$R_tilde, true_R = R_tilde)
  
  EV_42 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_42_tilde$L_tilde,
    G_tilde = kcv_42_tilde$G_tilde,
    R_tilde = kcv_42_tilde$R_tilde
  )
  
  ## 10-fold cross-validation--62
  kcv_62_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_62, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_62_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_62_opt$opt_para[1:2]), lambda=as.numeric(kcv_62_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_62_tilde <- MakeIdent(L=kcv_62_res$L, G=kcv_62_res$G, R=kcv_62_res$R)
  
  kcv_62_lossM <- loss(tnsr = kcv_62_res$est, true_tnsr=Msmooth@data)
  
  kcv_62_truncated_res <- rank_truncation(Mmiss, mode="L", true_rank_L=3, true_rank_R=NULL,
                                          L_tilde=kcv_62_tilde$L_tilde, G_tilde=kcv_62_tilde$G_tilde, R_tilde=kcv_62_tilde$R_tilde)
  
  kcv_62_loss_truncated <- loss(tnsr = kcv_62_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_62_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_62_tilde$R_tilde, true_R = R_tilde)
  
  EV_62 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_62_tilde$L_tilde,
    G_tilde = kcv_62_tilde$G_tilde,
    R_tilde = kcv_62_tilde$R_tilde
  )
  
  ## 10-fold cross-validation--33
  kcv_33_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_33, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_33_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_33_opt$opt_para[1:2]), lambda=as.numeric(kcv_33_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_33_tilde <- MakeIdent(L=kcv_33_res$L, G=kcv_33_res$G, R=kcv_33_res$R)
  
  kcv_33_lossM <- loss(tnsr = kcv_33_res$est, true_tnsr=Msmooth@data)
  
  kcv_33_truncated_res <- rank_truncation(Mmiss, mode="R", true_rank_L=NULL, true_rank_R=2,
                                          L_tilde=kcv_33_tilde$L_tilde, G_tilde=kcv_33_tilde$G_tilde, R_tilde=kcv_33_tilde$R_tilde)
  
  kcv_33_loss_truncated <- loss(tnsr = kcv_33_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_33_tilde$L_tilde, true_L = L_tilde,
                                R = kcv_33_truncated_res$R_truncated, true_R = R_tilde)
  
  EV_33 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_33_tilde$L_tilde,
    G_tilde = kcv_33_tilde$G_tilde,
    R_tilde = kcv_33_tilde$R_tilde
  )
  
  ## 10-fold cross-validation--43
  kcv_43_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_43, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_43_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_43_opt$opt_para[1:2]), lambda=as.numeric(kcv_43_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_43_tilde <- MakeIdent(L=kcv_43_res$L, G=kcv_43_res$G, R=kcv_43_res$R)
  
  kcv_43_lossM <- loss(tnsr = kcv_43_res$est, true_tnsr=Msmooth@data)
  
  kcv_43_truncated_res <- rank_truncation(Mmiss, mode="both", true_rank_L=3, true_rank_R=2,
                                          L_tilde=kcv_43_tilde$L_tilde, G_tilde=kcv_43_tilde$G_tilde, R_tilde=kcv_43_tilde$R_tilde)
  
  kcv_43_loss_truncated <- loss(tnsr = kcv_43_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_43_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_43_truncated_res$R_truncated, true_R = R_tilde)
  
  EV_43 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_43_tilde$L_tilde,
    G_tilde = kcv_43_tilde$G_tilde,
    R_tilde = kcv_43_tilde$R_tilde
  )
  
  ## 10-fold cross-validation--63
  kcv_63_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_63, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_63_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_63_opt$opt_para[1:2]), lambda=as.numeric(kcv_63_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_63_tilde <- MakeIdent(L=kcv_63_res$L, G=kcv_63_res$G, R=kcv_63_res$R)
  
  kcv_63_lossM <- loss(tnsr = kcv_63_res$est, true_tnsr=Msmooth@data)
  
  kcv_63_truncated_res <- rank_truncation(Mmiss, mode="both", true_rank_L=3, true_rank_R=2,
                                          L_tilde=kcv_63_tilde$L_tilde, G_tilde=kcv_63_tilde$G_tilde, R_tilde=kcv_63_tilde$R_tilde)
  
  kcv_63_loss_truncated <- loss(tnsr = kcv_63_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_63_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_63_truncated_res$R_truncated, true_R = R_tilde)
  
  
  EV_63 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_63_tilde$L_tilde,
    G_tilde = kcv_63_tilde$G_tilde,
    R_tilde = kcv_63_tilde$R_tilde
  )
  
 
  output <- list(Msim = Mmiss,
                 "kcv_true_loss"=kcv_true_loss, lambda_true = as.numeric(kcv_true_opt$opt_para[3]),
                 "kcv_42_lossM"=kcv_42_lossM, "kcv_42_loss_truncated"=kcv_42_loss_truncated,"lambda_42" = as.numeric(kcv_42_opt$opt_para[3]),
                 "kcv_62_lossM"=kcv_62_lossM, "kcv_62_loss_truncated"=kcv_62_loss_truncated,"lambda_62" = as.numeric(kcv_62_opt$opt_para[3]),
                 "kcv_33_lossM"=kcv_33_lossM, "kcv_33_loss_truncated"=kcv_33_loss_truncated,"lambda_33" = as.numeric(kcv_33_opt$opt_para[3]),
                 "kcv_43_lossM"=kcv_43_lossM, "kcv_43_loss_truncated"=kcv_43_loss_truncated,"lambda_43" = as.numeric(kcv_43_opt$opt_para[3]),
                 "kcv_63_lossM"=kcv_63_lossM, "kcv_63_loss_truncated"=kcv_63_loss_truncated,"lambda_63" = as.numeric(kcv_63_opt$opt_para[3]),
                 "EV_true" = EV_true,"EV_42" = EV_42,"EV_62" = EV_62, "EV_33" = EV_33,"EV_43" = EV_43,"EV_63" = EV_63
  )
  output
}

save(miss20, file="./miss20_ev.Rda")

## missing rate 50%
set.seed(20261004)

miss50 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","multiway","funData","MFPCA","SmoothHOOI")) %dorng%{
  sim_data <- sim_data2(L = L_tilde, R= R_tilde, mean_G = mean_G, cov_G = cov_G,
                        E = E, p=200, noise_level=1, pattern="random", percent=0.5)
  Mmiss <- sim_data$sim_Mmiss
  Msmooth <- sim_data$sim_Msmooth
  
  ## 10-fold cross-validation--true rank
  kcv_true_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=true_rank_grid, lambda_seq=lambda_seq, k=10,
                             L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_true_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_true_opt$opt_para[1:2]), lambda=as.numeric(kcv_true_opt$opt_para[3]),
                         L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_true_tilde <- MakeIdent(L=kcv_true_res$L, G=kcv_true_res$G, R=kcv_true_res$R)
  
  kcv_true_loss <- loss(tnsr = kcv_true_res$est, true_tnsr=Msmooth@data,
                        L = kcv_true_tilde$L_tilde, true_L = L_tilde,
                        R = kcv_true_tilde$R_tilde, true_R = R_tilde)
  
  
  EV_true <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_true_tilde$L_tilde,
    G_tilde = kcv_true_tilde$G_tilde,
    R_tilde = kcv_true_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--42
  kcv_42_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_42, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_42_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_42_opt$opt_para[1:2]), lambda=as.numeric(kcv_42_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_42_tilde <- MakeIdent(L=kcv_42_res$L, G=kcv_42_res$G, R=kcv_42_res$R)
  
  kcv_42_lossM <- loss(tnsr = kcv_42_res$est, true_tnsr=Msmooth@data)
  
  kcv_42_truncated_res <- rank_truncation(Mmiss, mode="L", true_rank_L=3, true_rank_R=NULL,
                                          L_tilde=kcv_42_tilde$L_tilde, G_tilde=kcv_42_tilde$G_tilde, R_tilde=kcv_42_tilde$R_tilde)
  
  kcv_42_loss_truncated <- loss(tnsr = kcv_42_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_42_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_42_tilde$R_tilde, true_R = R_tilde)
  
  EV_42 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_42_tilde$L_tilde,
    G_tilde = kcv_42_tilde$G_tilde,
    R_tilde = kcv_42_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--62
  kcv_62_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_62, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_62_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_62_opt$opt_para[1:2]), lambda=as.numeric(kcv_62_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_62_tilde <- MakeIdent(L=kcv_62_res$L, G=kcv_62_res$G, R=kcv_62_res$R)
  
  kcv_62_lossM <- loss(tnsr = kcv_62_res$est, true_tnsr=Msmooth@data)
  
  kcv_62_truncated_res <- rank_truncation(Mmiss, mode="L", true_rank_L=3, true_rank_R=NULL,
                                          L_tilde=kcv_62_tilde$L_tilde, G_tilde=kcv_62_tilde$G_tilde, R_tilde=kcv_62_tilde$R_tilde)
  
  kcv_62_loss_truncated <- loss(tnsr = kcv_62_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_62_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_62_tilde$R_tilde, true_R = R_tilde)
  
  EV_62 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_62_tilde$L_tilde,
    G_tilde = kcv_62_tilde$G_tilde,
    R_tilde = kcv_62_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--33
  kcv_33_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_33, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_33_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_33_opt$opt_para[1:2]), lambda=as.numeric(kcv_33_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_33_tilde <- MakeIdent(L=kcv_33_res$L, G=kcv_33_res$G, R=kcv_33_res$R)
  
  kcv_33_lossM <- loss(tnsr = kcv_33_res$est, true_tnsr=Msmooth@data)
  
  kcv_33_truncated_res <- rank_truncation(Mmiss, mode="R", true_rank_L=NULL, true_rank_R=2,
                                          L_tilde=kcv_33_tilde$L_tilde, G_tilde=kcv_33_tilde$G_tilde, R_tilde=kcv_33_tilde$R_tilde)
  
  kcv_33_loss_truncated <- loss(tnsr = kcv_33_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_33_tilde$L_tilde, true_L = L_tilde,
                                R = kcv_33_truncated_res$R_truncated, true_R = R_tilde)
  
  EV_33 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_33_tilde$L_tilde,
    G_tilde = kcv_33_tilde$G_tilde,
    R_tilde = kcv_33_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--43
  kcv_43_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_43, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_43_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_43_opt$opt_para[1:2]), lambda=as.numeric(kcv_43_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_43_tilde <- MakeIdent(L=kcv_43_res$L, G=kcv_43_res$G, R=kcv_43_res$R)
  
  kcv_43_lossM <- loss(tnsr = kcv_43_res$est, true_tnsr=Msmooth@data)
  
  kcv_43_truncated_res <- rank_truncation(Mmiss, mode="both", true_rank_L=3, true_rank_R=2,
                                          L_tilde=kcv_43_tilde$L_tilde, G_tilde=kcv_43_tilde$G_tilde, R_tilde=kcv_43_tilde$R_tilde)
  
  kcv_43_loss_truncated <- loss(tnsr = kcv_43_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_43_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_43_truncated_res$R_truncated, true_R = R_tilde)
  
  EV_43 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_43_tilde$L_tilde,
    G_tilde = kcv_43_tilde$G_tilde,
    R_tilde = kcv_43_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--63
  kcv_63_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_63, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_63_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_63_opt$opt_para[1:2]), lambda=as.numeric(kcv_63_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_63_tilde <- MakeIdent(L=kcv_63_res$L, G=kcv_63_res$G, R=kcv_63_res$R)
  
  kcv_63_lossM <- loss(tnsr = kcv_63_res$est, true_tnsr=Msmooth@data)
  
  kcv_63_truncated_res <- rank_truncation(Mmiss, mode="both", true_rank_L=3, true_rank_R=2,
                                          L_tilde=kcv_63_tilde$L_tilde, G_tilde=kcv_63_tilde$G_tilde, R_tilde=kcv_63_tilde$R_tilde)
  
  kcv_63_loss_truncated <- loss(tnsr = kcv_63_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_63_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_63_truncated_res$R_truncated, true_R = R_tilde)
  
  
  EV_63 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_63_tilde$L_tilde,
    G_tilde = kcv_63_tilde$G_tilde,
    R_tilde = kcv_63_tilde$R_tilde
  )
  
  output <- list("kcv_true_loss"=kcv_true_loss, lambda_true = as.numeric(kcv_true_opt$opt_para[3]),
                 "kcv_42_lossM"=kcv_42_lossM, "kcv_42_loss_truncated"=kcv_42_loss_truncated,"lambda_42" = as.numeric(kcv_42_opt$opt_para[3]),
                 "kcv_62_lossM"=kcv_62_lossM, "kcv_62_loss_truncated"=kcv_62_loss_truncated,"lambda_62" = as.numeric(kcv_62_opt$opt_para[3]),
                 "kcv_33_lossM"=kcv_33_lossM, "kcv_33_loss_truncated"=kcv_33_loss_truncated,"lambda_33" = as.numeric(kcv_33_opt$opt_para[3]),
                 "kcv_43_lossM"=kcv_43_lossM, "kcv_43_loss_truncated"=kcv_43_loss_truncated,"lambda_43" = as.numeric(kcv_43_opt$opt_para[3]),
                 "kcv_63_lossM"=kcv_63_lossM, "kcv_63_loss_truncated"=kcv_63_loss_truncated,"lambda_63" = as.numeric(kcv_63_opt$opt_para[3]),
                 "EV_true" = EV_true,"EV_42" = EV_42,"EV_62" = EV_62, "EV_33" = EV_33,"EV_43" = EV_43,"EV_63" = EV_63
  )
  output
}

save(miss50, file="./rank_misspecification_ev/miss50.Rda")

## structured missing
set.seed(20261005)

missstruc <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","multiway","funData","MFPCA","SmoothHOOI")) %dorng%{
  sim_data <- sim_data2(L = L_tilde, R= R_tilde, mean_G = mean_G, cov_G = cov_G,
                        E = E, p=200, noise_level=1, pattern="structured", lower=0, upper=20)
  Mmiss <- sim_data$sim_Mmiss
  Msmooth <- sim_data$sim_Msmooth
  
  missing_rate <- sum(is.na(Mmiss@data))/(prod(dim(Mmiss@data)))
  
  ## 10-fold cross-validation--true rank
  kcv_true_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=true_rank_grid, lambda_seq=lambda_seq, k=10,
                             L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_true_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_true_opt$opt_para[1:2]), lambda=as.numeric(kcv_true_opt$opt_para[3]),
                         L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_true_tilde <- MakeIdent(L=kcv_true_res$L, G=kcv_true_res$G, R=kcv_true_res$R)
  
  kcv_true_loss <- loss(tnsr = kcv_true_res$est, true_tnsr=Msmooth@data,
                        L = kcv_true_tilde$L_tilde, true_L = L_tilde,
                        R = kcv_true_tilde$R_tilde, true_R = R_tilde)
  
  
  EV_true <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_true_tilde$L_tilde,
    G_tilde = kcv_true_tilde$G_tilde,
    R_tilde = kcv_true_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--42
  kcv_42_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_42, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_42_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_42_opt$opt_para[1:2]), lambda=as.numeric(kcv_42_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_42_tilde <- MakeIdent(L=kcv_42_res$L, G=kcv_42_res$G, R=kcv_42_res$R)
  
  kcv_42_lossM <- loss(tnsr = kcv_42_res$est, true_tnsr=Msmooth@data)
  
  kcv_42_truncated_res <- rank_truncation(Mmiss, mode="L", true_rank_L=3, true_rank_R=NULL,
                                          L_tilde=kcv_42_tilde$L_tilde, G_tilde=kcv_42_tilde$G_tilde, R_tilde=kcv_42_tilde$R_tilde)
  
  kcv_42_loss_truncated <- loss(tnsr = kcv_42_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_42_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_42_tilde$R_tilde, true_R = R_tilde)
  
  EV_42 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_42_tilde$L_tilde,
    G_tilde = kcv_42_tilde$G_tilde,
    R_tilde = kcv_42_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--62
  kcv_62_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_62, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_62_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_62_opt$opt_para[1:2]), lambda=as.numeric(kcv_62_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_62_tilde <- MakeIdent(L=kcv_62_res$L, G=kcv_62_res$G, R=kcv_62_res$R)
  
  kcv_62_lossM <- loss(tnsr = kcv_62_res$est, true_tnsr=Msmooth@data)
  
  kcv_62_truncated_res <- rank_truncation(Mmiss, mode="L", true_rank_L=3, true_rank_R=NULL,
                                          L_tilde=kcv_62_tilde$L_tilde, G_tilde=kcv_62_tilde$G_tilde, R_tilde=kcv_62_tilde$R_tilde)
  
  kcv_62_loss_truncated <- loss(tnsr = kcv_62_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_62_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_62_tilde$R_tilde, true_R = R_tilde)
  
  EV_62 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_62_tilde$L_tilde,
    G_tilde = kcv_62_tilde$G_tilde,
    R_tilde = kcv_62_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--33
  kcv_33_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_33, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_33_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_33_opt$opt_para[1:2]), lambda=as.numeric(kcv_33_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_33_tilde <- MakeIdent(L=kcv_33_res$L, G=kcv_33_res$G, R=kcv_33_res$R)
  
  kcv_33_lossM <- loss(tnsr = kcv_33_res$est, true_tnsr=Msmooth@data)
  
  kcv_33_truncated_res <- rank_truncation(Mmiss, mode="R", true_rank_L=NULL, true_rank_R=2,
                                          L_tilde=kcv_33_tilde$L_tilde, G_tilde=kcv_33_tilde$G_tilde, R_tilde=kcv_33_tilde$R_tilde)
  
  kcv_33_loss_truncated <- loss(tnsr = kcv_33_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_33_tilde$L_tilde, true_L = L_tilde,
                                R = kcv_33_truncated_res$R_truncated, true_R = R_tilde)
  
  EV_33 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_33_tilde$L_tilde,
    G_tilde = kcv_33_tilde$G_tilde,
    R_tilde = kcv_33_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--43
  kcv_43_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_43, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_43_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_43_opt$opt_para[1:2]), lambda=as.numeric(kcv_43_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_43_tilde <- MakeIdent(L=kcv_43_res$L, G=kcv_43_res$G, R=kcv_43_res$R)
  
  kcv_43_lossM <- loss(tnsr = kcv_43_res$est, true_tnsr=Msmooth@data)
  
  kcv_43_truncated_res <- rank_truncation(Mmiss, mode="both", true_rank_L=3, true_rank_R=2,
                                          L_tilde=kcv_43_tilde$L_tilde, G_tilde=kcv_43_tilde$G_tilde, R_tilde=kcv_43_tilde$R_tilde)
  
  kcv_43_loss_truncated <- loss(tnsr = kcv_43_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_43_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_43_truncated_res$R_truncated, true_R = R_tilde)
  
  EV_43 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_43_tilde$L_tilde,
    G_tilde = kcv_43_tilde$G_tilde,
    R_tilde = kcv_43_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--63
  kcv_63_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_63, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_63_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_63_opt$opt_para[1:2]), lambda=as.numeric(kcv_63_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_63_tilde <- MakeIdent(L=kcv_63_res$L, G=kcv_63_res$G, R=kcv_63_res$R)
  
  kcv_63_lossM <- loss(tnsr = kcv_63_res$est, true_tnsr=Msmooth@data)
  
  kcv_63_truncated_res <- rank_truncation(Mmiss, mode="both", true_rank_L=3, true_rank_R=2,
                                          L_tilde=kcv_63_tilde$L_tilde, G_tilde=kcv_63_tilde$G_tilde, R_tilde=kcv_63_tilde$R_tilde)
  
  kcv_63_loss_truncated <- loss(tnsr = kcv_63_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_63_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_63_truncated_res$R_truncated, true_R = R_tilde)
  
  
  EV_63 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_63_tilde$L_tilde,
    G_tilde = kcv_63_tilde$G_tilde,
    R_tilde = kcv_63_tilde$R_tilde
  )
  
  output <- list("kcv_true_loss"=kcv_true_loss, lambda_true = as.numeric(kcv_true_opt$opt_para[3]),
                 "kcv_42_lossM"=kcv_42_lossM, "kcv_42_loss_truncated"=kcv_42_loss_truncated,"lambda_42" = as.numeric(kcv_42_opt$opt_para[3]),
                 "kcv_62_lossM"=kcv_62_lossM, "kcv_62_loss_truncated"=kcv_62_loss_truncated,"lambda_62" = as.numeric(kcv_62_opt$opt_para[3]),
                 "kcv_33_lossM"=kcv_33_lossM, "kcv_33_loss_truncated"=kcv_33_loss_truncated,"lambda_33" = as.numeric(kcv_33_opt$opt_para[3]),
                 "kcv_43_lossM"=kcv_43_lossM, "kcv_43_loss_truncated"=kcv_43_loss_truncated,"lambda_43" = as.numeric(kcv_43_opt$opt_para[3]),
                 "kcv_63_lossM"=kcv_63_lossM, "kcv_63_loss_truncated"=kcv_63_loss_truncated,"lambda_63" = as.numeric(kcv_63_opt$opt_para[3]),
                 "EV_true" = EV_true,"EV_42" = EV_42,"EV_62" = EV_62, "EV_33" = EV_33,"EV_43" = EV_43,"EV_63" = EV_63,
                 missing_rate=missing_rate
  )
  output
}

save(missstruc, file="./rank_misspecification_ev/missstruc.Rda")


## p = 50
set.seed(20261006)

p50 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","multiway","funData","MFPCA","SmoothHOOI")) %dorng%{
  sim_data <- sim_data2(L = L_tilde, R= R_tilde, mean_G = mean_G, cov_G = cov_G,
                        E = E, p=50, noise_level=1, pattern="random", percent=0.2)
  Mmiss <- sim_data$sim_Mmiss
  Msmooth <- sim_data$sim_Msmooth
  
  ## 10-fold cross-validation--true rank
  kcv_true_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=true_rank_grid, lambda_seq=lambda_seq, k=10,
                             L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_true_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_true_opt$opt_para[1:2]), lambda=as.numeric(kcv_true_opt$opt_para[3]),
                         L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_true_tilde <- MakeIdent(L=kcv_true_res$L, G=kcv_true_res$G, R=kcv_true_res$R)
  
  kcv_true_loss <- loss(tnsr = kcv_true_res$est, true_tnsr=Msmooth@data,
                        L = kcv_true_tilde$L_tilde, true_L = L_tilde,
                        R = kcv_true_tilde$R_tilde, true_R = R_tilde)
  
  
  EV_true <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_true_tilde$L_tilde,
    G_tilde = kcv_true_tilde$G_tilde,
    R_tilde = kcv_true_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--42
  kcv_42_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_42, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_42_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_42_opt$opt_para[1:2]), lambda=as.numeric(kcv_42_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_42_tilde <- MakeIdent(L=kcv_42_res$L, G=kcv_42_res$G, R=kcv_42_res$R)
  
  kcv_42_lossM <- loss(tnsr = kcv_42_res$est, true_tnsr=Msmooth@data)
  
  kcv_42_truncated_res <- rank_truncation(Mmiss, mode="L", true_rank_L=3, true_rank_R=NULL,
                                          L_tilde=kcv_42_tilde$L_tilde, G_tilde=kcv_42_tilde$G_tilde, R_tilde=kcv_42_tilde$R_tilde)
  
  kcv_42_loss_truncated <- loss(tnsr = kcv_42_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_42_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_42_tilde$R_tilde, true_R = R_tilde)
  
  EV_42 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_42_tilde$L_tilde,
    G_tilde = kcv_42_tilde$G_tilde,
    R_tilde = kcv_42_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--62
  kcv_62_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_62, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_62_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_62_opt$opt_para[1:2]), lambda=as.numeric(kcv_62_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_62_tilde <- MakeIdent(L=kcv_62_res$L, G=kcv_62_res$G, R=kcv_62_res$R)
  
  kcv_62_lossM <- loss(tnsr = kcv_62_res$est, true_tnsr=Msmooth@data)
  
  kcv_62_truncated_res <- rank_truncation(Mmiss, mode="L", true_rank_L=3, true_rank_R=NULL,
                                          L_tilde=kcv_62_tilde$L_tilde, G_tilde=kcv_62_tilde$G_tilde, R_tilde=kcv_62_tilde$R_tilde)
  
  kcv_62_loss_truncated <- loss(tnsr = kcv_62_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_62_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_62_tilde$R_tilde, true_R = R_tilde)
  
  EV_62 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_62_tilde$L_tilde,
    G_tilde = kcv_62_tilde$G_tilde,
    R_tilde = kcv_62_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--33
  kcv_33_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_33, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_33_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_33_opt$opt_para[1:2]), lambda=as.numeric(kcv_33_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_33_tilde <- MakeIdent(L=kcv_33_res$L, G=kcv_33_res$G, R=kcv_33_res$R)
  
  kcv_33_lossM <- loss(tnsr = kcv_33_res$est, true_tnsr=Msmooth@data)
  
  kcv_33_truncated_res <- rank_truncation(Mmiss, mode="R", true_rank_L=NULL, true_rank_R=2,
                                          L_tilde=kcv_33_tilde$L_tilde, G_tilde=kcv_33_tilde$G_tilde, R_tilde=kcv_33_tilde$R_tilde)
  
  kcv_33_loss_truncated <- loss(tnsr = kcv_33_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_33_tilde$L_tilde, true_L = L_tilde,
                                R = kcv_33_truncated_res$R_truncated, true_R = R_tilde)
  
  EV_33 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_33_tilde$L_tilde,
    G_tilde = kcv_33_tilde$G_tilde,
    R_tilde = kcv_33_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--43
  kcv_43_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_43, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_43_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_43_opt$opt_para[1:2]), lambda=as.numeric(kcv_43_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_43_tilde <- MakeIdent(L=kcv_43_res$L, G=kcv_43_res$G, R=kcv_43_res$R)
  
  kcv_43_lossM <- loss(tnsr = kcv_43_res$est, true_tnsr=Msmooth@data)
  
  kcv_43_truncated_res <- rank_truncation(Mmiss, mode="both", true_rank_L=3, true_rank_R=2,
                                          L_tilde=kcv_43_tilde$L_tilde, G_tilde=kcv_43_tilde$G_tilde, R_tilde=kcv_43_tilde$R_tilde)
  
  kcv_43_loss_truncated <- loss(tnsr = kcv_43_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_43_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_43_truncated_res$R_truncated, true_R = R_tilde)
  
  EV_43 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_43_tilde$L_tilde,
    G_tilde = kcv_43_tilde$G_tilde,
    R_tilde = kcv_43_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--63
  kcv_63_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_63, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_63_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_63_opt$opt_para[1:2]), lambda=as.numeric(kcv_63_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_63_tilde <- MakeIdent(L=kcv_63_res$L, G=kcv_63_res$G, R=kcv_63_res$R)
  
  kcv_63_lossM <- loss(tnsr = kcv_63_res$est, true_tnsr=Msmooth@data)
  
  kcv_63_truncated_res <- rank_truncation(Mmiss, mode="both", true_rank_L=3, true_rank_R=2,
                                          L_tilde=kcv_63_tilde$L_tilde, G_tilde=kcv_63_tilde$G_tilde, R_tilde=kcv_63_tilde$R_tilde)
  
  kcv_63_loss_truncated <- loss(tnsr = kcv_63_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_63_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_63_truncated_res$R_truncated, true_R = R_tilde)
  
  
  EV_63 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_63_tilde$L_tilde,
    G_tilde = kcv_63_tilde$G_tilde,
    R_tilde = kcv_63_tilde$R_tilde
  )
  
  output <- list("kcv_true_loss"=kcv_true_loss, lambda_true = as.numeric(kcv_true_opt$opt_para[3]),
                 "kcv_42_lossM"=kcv_42_lossM, "kcv_42_loss_truncated"=kcv_42_loss_truncated,"lambda_42" = as.numeric(kcv_42_opt$opt_para[3]),
                 "kcv_62_lossM"=kcv_62_lossM, "kcv_62_loss_truncated"=kcv_62_loss_truncated,"lambda_62" = as.numeric(kcv_62_opt$opt_para[3]),
                 "kcv_33_lossM"=kcv_33_lossM, "kcv_33_loss_truncated"=kcv_33_loss_truncated,"lambda_33" = as.numeric(kcv_33_opt$opt_para[3]),
                 "kcv_43_lossM"=kcv_43_lossM, "kcv_43_loss_truncated"=kcv_43_loss_truncated,"lambda_43" = as.numeric(kcv_43_opt$opt_para[3]),
                 "kcv_63_lossM"=kcv_63_lossM, "kcv_63_loss_truncated"=kcv_63_loss_truncated,"lambda_63" = as.numeric(kcv_63_opt$opt_para[3]),
                 "EV_true" = EV_true,"EV_42" = EV_42,"EV_62" = EV_62, "EV_33" = EV_33,"EV_43" = EV_43,"EV_63" = EV_63
  )
  output
}

save(p50, file="./rank_misspecification_ev/p50.Rda")

## p = 500
set.seed(20261007)

p500 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","multiway","funData","MFPCA","SmoothHOOI")) %dorng%{
  sim_data <- sim_data2(L = L_tilde, R= R_tilde, mean_G = mean_G, cov_G = cov_G,
                        E = E, p=500, noise_level=1, pattern="random", percent=0.2)
  Mmiss <- sim_data$sim_Mmiss
  Msmooth <- sim_data$sim_Msmooth
  
  ## 10-fold cross-validation--true rank
  kcv_true_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=true_rank_grid, lambda_seq=lambda_seq, k=10,
                             L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_true_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_true_opt$opt_para[1:2]), lambda=as.numeric(kcv_true_opt$opt_para[3]),
                         L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_true_tilde <- MakeIdent(L=kcv_true_res$L, G=kcv_true_res$G, R=kcv_true_res$R)
  
  kcv_true_loss <- loss(tnsr = kcv_true_res$est, true_tnsr=Msmooth@data,
                        L = kcv_true_tilde$L_tilde, true_L = L_tilde,
                        R = kcv_true_tilde$R_tilde, true_R = R_tilde)
  
  
  EV_true <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_true_tilde$L_tilde,
    G_tilde = kcv_true_tilde$G_tilde,
    R_tilde = kcv_true_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--42
  kcv_42_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_42, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_42_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_42_opt$opt_para[1:2]), lambda=as.numeric(kcv_42_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_42_tilde <- MakeIdent(L=kcv_42_res$L, G=kcv_42_res$G, R=kcv_42_res$R)
  
  kcv_42_lossM <- loss(tnsr = kcv_42_res$est, true_tnsr=Msmooth@data)
  
  kcv_42_truncated_res <- rank_truncation(Mmiss, mode="L", true_rank_L=3, true_rank_R=NULL,
                                          L_tilde=kcv_42_tilde$L_tilde, G_tilde=kcv_42_tilde$G_tilde, R_tilde=kcv_42_tilde$R_tilde)
  
  kcv_42_loss_truncated <- loss(tnsr = kcv_42_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_42_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_42_tilde$R_tilde, true_R = R_tilde)
  
  EV_42 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_42_tilde$L_tilde,
    G_tilde = kcv_42_tilde$G_tilde,
    R_tilde = kcv_42_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--62
  kcv_62_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_62, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_62_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_62_opt$opt_para[1:2]), lambda=as.numeric(kcv_62_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_62_tilde <- MakeIdent(L=kcv_62_res$L, G=kcv_62_res$G, R=kcv_62_res$R)
  
  kcv_62_lossM <- loss(tnsr = kcv_62_res$est, true_tnsr=Msmooth@data)
  
  kcv_62_truncated_res <- rank_truncation(Mmiss, mode="L", true_rank_L=3, true_rank_R=NULL,
                                          L_tilde=kcv_62_tilde$L_tilde, G_tilde=kcv_62_tilde$G_tilde, R_tilde=kcv_62_tilde$R_tilde)
  
  kcv_62_loss_truncated <- loss(tnsr = kcv_62_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_62_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_62_tilde$R_tilde, true_R = R_tilde)
  
  EV_62 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_62_tilde$L_tilde,
    G_tilde = kcv_62_tilde$G_tilde,
    R_tilde = kcv_62_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--33
  kcv_33_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_33, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_33_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_33_opt$opt_para[1:2]), lambda=as.numeric(kcv_33_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_33_tilde <- MakeIdent(L=kcv_33_res$L, G=kcv_33_res$G, R=kcv_33_res$R)
  
  kcv_33_lossM <- loss(tnsr = kcv_33_res$est, true_tnsr=Msmooth@data)
  
  kcv_33_truncated_res <- rank_truncation(Mmiss, mode="R", true_rank_L=NULL, true_rank_R=2,
                                          L_tilde=kcv_33_tilde$L_tilde, G_tilde=kcv_33_tilde$G_tilde, R_tilde=kcv_33_tilde$R_tilde)
  
  kcv_33_loss_truncated <- loss(tnsr = kcv_33_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_33_tilde$L_tilde, true_L = L_tilde,
                                R = kcv_33_truncated_res$R_truncated, true_R = R_tilde)
  
  EV_33 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_33_tilde$L_tilde,
    G_tilde = kcv_33_tilde$G_tilde,
    R_tilde = kcv_33_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--43
  kcv_43_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_43, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_43_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_43_opt$opt_para[1:2]), lambda=as.numeric(kcv_43_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_43_tilde <- MakeIdent(L=kcv_43_res$L, G=kcv_43_res$G, R=kcv_43_res$R)
  
  kcv_43_lossM <- loss(tnsr = kcv_43_res$est, true_tnsr=Msmooth@data)
  
  kcv_43_truncated_res <- rank_truncation(Mmiss, mode="both", true_rank_L=3, true_rank_R=2,
                                          L_tilde=kcv_43_tilde$L_tilde, G_tilde=kcv_43_tilde$G_tilde, R_tilde=kcv_43_tilde$R_tilde)
  
  kcv_43_loss_truncated <- loss(tnsr = kcv_43_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_43_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_43_truncated_res$R_truncated, true_R = R_tilde)
  
  EV_43 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_43_tilde$L_tilde,
    G_tilde = kcv_43_tilde$G_tilde,
    R_tilde = kcv_43_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--63
  kcv_63_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_63, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_63_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_63_opt$opt_para[1:2]), lambda=as.numeric(kcv_63_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_63_tilde <- MakeIdent(L=kcv_63_res$L, G=kcv_63_res$G, R=kcv_63_res$R)
  
  kcv_63_lossM <- loss(tnsr = kcv_63_res$est, true_tnsr=Msmooth@data)
  
  kcv_63_truncated_res <- rank_truncation(Mmiss, mode="both", true_rank_L=3, true_rank_R=2,
                                          L_tilde=kcv_63_tilde$L_tilde, G_tilde=kcv_63_tilde$G_tilde, R_tilde=kcv_63_tilde$R_tilde)
  
  kcv_63_loss_truncated <- loss(tnsr = kcv_63_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_63_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_63_truncated_res$R_truncated, true_R = R_tilde)
  
  
  EV_63 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_63_tilde$L_tilde,
    G_tilde = kcv_63_tilde$G_tilde,
    R_tilde = kcv_63_tilde$R_tilde
  )
  
  output <- list("kcv_true_loss"=kcv_true_loss, lambda_true = as.numeric(kcv_true_opt$opt_para[3]),
                 "kcv_42_lossM"=kcv_42_lossM, "kcv_42_loss_truncated"=kcv_42_loss_truncated,"lambda_42" = as.numeric(kcv_42_opt$opt_para[3]),
                 "kcv_62_lossM"=kcv_62_lossM, "kcv_62_loss_truncated"=kcv_62_loss_truncated,"lambda_62" = as.numeric(kcv_62_opt$opt_para[3]),
                 "kcv_33_lossM"=kcv_33_lossM, "kcv_33_loss_truncated"=kcv_33_loss_truncated,"lambda_33" = as.numeric(kcv_33_opt$opt_para[3]),
                 "kcv_43_lossM"=kcv_43_lossM, "kcv_43_loss_truncated"=kcv_43_loss_truncated,"lambda_43" = as.numeric(kcv_43_opt$opt_para[3]),
                 "kcv_63_lossM"=kcv_63_lossM, "kcv_63_loss_truncated"=kcv_63_loss_truncated,"lambda_63" = as.numeric(kcv_63_opt$opt_para[3]),
                 "EV_true" = EV_true,"EV_42" = EV_42,"EV_62" = EV_62, "EV_33" = EV_33,"EV_43" = EV_43,"EV_63" = EV_63
  )
  output
}

save(p500, file="./rank_misspecification_ev/p500.Rda")


## noise level 0*empirical var
set.seed(20261008)

noise0 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","multiway","funData","MFPCA","SmoothHOOI")) %dorng%{
  sim_data <- sim_data2(L = L_tilde, R= R_tilde, mean_G = mean_G, cov_G = cov_G,
                        E = E, p=200, noise_level=0, pattern="random", percent=0.2)
  Mmiss <- sim_data$sim_Mmiss
  Msmooth <- sim_data$sim_Msmooth
  
  ## 10-fold cross-validation--true rank
  kcv_true_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=true_rank_grid, lambda_seq=lambda_seq, k=10,
                             L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_true_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_true_opt$opt_para[1:2]), lambda=as.numeric(kcv_true_opt$opt_para[3]),
                         L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_true_tilde <- MakeIdent(L=kcv_true_res$L, G=kcv_true_res$G, R=kcv_true_res$R)
  
  kcv_true_loss <- loss(tnsr = kcv_true_res$est, true_tnsr=Msmooth@data,
                        L = kcv_true_tilde$L_tilde, true_L = L_tilde,
                        R = kcv_true_tilde$R_tilde, true_R = R_tilde)
  
  
  EV_true <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_true_tilde$L_tilde,
    G_tilde = kcv_true_tilde$G_tilde,
    R_tilde = kcv_true_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--42
  kcv_42_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_42, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_42_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_42_opt$opt_para[1:2]), lambda=as.numeric(kcv_42_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_42_tilde <- MakeIdent(L=kcv_42_res$L, G=kcv_42_res$G, R=kcv_42_res$R)
  
  kcv_42_lossM <- loss(tnsr = kcv_42_res$est, true_tnsr=Msmooth@data)
  
  kcv_42_truncated_res <- rank_truncation(Mmiss, mode="L", true_rank_L=3, true_rank_R=NULL,
                                          L_tilde=kcv_42_tilde$L_tilde, G_tilde=kcv_42_tilde$G_tilde, R_tilde=kcv_42_tilde$R_tilde)
  
  kcv_42_loss_truncated <- loss(tnsr = kcv_42_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_42_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_42_tilde$R_tilde, true_R = R_tilde)
  
  EV_42 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_42_tilde$L_tilde,
    G_tilde = kcv_42_tilde$G_tilde,
    R_tilde = kcv_42_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--62
  kcv_62_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_62, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_62_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_62_opt$opt_para[1:2]), lambda=as.numeric(kcv_62_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_62_tilde <- MakeIdent(L=kcv_62_res$L, G=kcv_62_res$G, R=kcv_62_res$R)
  
  kcv_62_lossM <- loss(tnsr = kcv_62_res$est, true_tnsr=Msmooth@data)
  
  kcv_62_truncated_res <- rank_truncation(Mmiss, mode="L", true_rank_L=3, true_rank_R=NULL,
                                          L_tilde=kcv_62_tilde$L_tilde, G_tilde=kcv_62_tilde$G_tilde, R_tilde=kcv_62_tilde$R_tilde)
  
  kcv_62_loss_truncated <- loss(tnsr = kcv_62_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_62_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_62_tilde$R_tilde, true_R = R_tilde)
  
  EV_62 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_62_tilde$L_tilde,
    G_tilde = kcv_62_tilde$G_tilde,
    R_tilde = kcv_62_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--33
  kcv_33_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_33, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_33_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_33_opt$opt_para[1:2]), lambda=as.numeric(kcv_33_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_33_tilde <- MakeIdent(L=kcv_33_res$L, G=kcv_33_res$G, R=kcv_33_res$R)
  
  kcv_33_lossM <- loss(tnsr = kcv_33_res$est, true_tnsr=Msmooth@data)
  
  kcv_33_truncated_res <- rank_truncation(Mmiss, mode="R", true_rank_L=NULL, true_rank_R=2,
                                          L_tilde=kcv_33_tilde$L_tilde, G_tilde=kcv_33_tilde$G_tilde, R_tilde=kcv_33_tilde$R_tilde)
  
  kcv_33_loss_truncated <- loss(tnsr = kcv_33_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_33_tilde$L_tilde, true_L = L_tilde,
                                R = kcv_33_truncated_res$R_truncated, true_R = R_tilde)
  
  EV_33 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_33_tilde$L_tilde,
    G_tilde = kcv_33_tilde$G_tilde,
    R_tilde = kcv_33_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--43
  kcv_43_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_43, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_43_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_43_opt$opt_para[1:2]), lambda=as.numeric(kcv_43_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_43_tilde <- MakeIdent(L=kcv_43_res$L, G=kcv_43_res$G, R=kcv_43_res$R)
  
  kcv_43_lossM <- loss(tnsr = kcv_43_res$est, true_tnsr=Msmooth@data)
  
  kcv_43_truncated_res <- rank_truncation(Mmiss, mode="both", true_rank_L=3, true_rank_R=2,
                                          L_tilde=kcv_43_tilde$L_tilde, G_tilde=kcv_43_tilde$G_tilde, R_tilde=kcv_43_tilde$R_tilde)
  
  kcv_43_loss_truncated <- loss(tnsr = kcv_43_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_43_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_43_truncated_res$R_truncated, true_R = R_tilde)
  
  EV_43 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_43_tilde$L_tilde,
    G_tilde = kcv_43_tilde$G_tilde,
    R_tilde = kcv_43_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--63
  kcv_63_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_63, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_63_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_63_opt$opt_para[1:2]), lambda=as.numeric(kcv_63_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_63_tilde <- MakeIdent(L=kcv_63_res$L, G=kcv_63_res$G, R=kcv_63_res$R)
  
  kcv_63_lossM <- loss(tnsr = kcv_63_res$est, true_tnsr=Msmooth@data)
  
  kcv_63_truncated_res <- rank_truncation(Mmiss, mode="both", true_rank_L=3, true_rank_R=2,
                                          L_tilde=kcv_63_tilde$L_tilde, G_tilde=kcv_63_tilde$G_tilde, R_tilde=kcv_63_tilde$R_tilde)
  
  kcv_63_loss_truncated <- loss(tnsr = kcv_63_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_63_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_63_truncated_res$R_truncated, true_R = R_tilde)
  
  
  EV_63 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_63_tilde$L_tilde,
    G_tilde = kcv_63_tilde$G_tilde,
    R_tilde = kcv_63_tilde$R_tilde
  )
  
  output <- list("kcv_true_loss"=kcv_true_loss, lambda_true = as.numeric(kcv_true_opt$opt_para[3]),
                 "kcv_42_lossM"=kcv_42_lossM, "kcv_42_loss_truncated"=kcv_42_loss_truncated,"lambda_42" = as.numeric(kcv_42_opt$opt_para[3]),
                 "kcv_62_lossM"=kcv_62_lossM, "kcv_62_loss_truncated"=kcv_62_loss_truncated,"lambda_62" = as.numeric(kcv_62_opt$opt_para[3]),
                 "kcv_33_lossM"=kcv_33_lossM, "kcv_33_loss_truncated"=kcv_33_loss_truncated,"lambda_33" = as.numeric(kcv_33_opt$opt_para[3]),
                 "kcv_43_lossM"=kcv_43_lossM, "kcv_43_loss_truncated"=kcv_43_loss_truncated,"lambda_43" = as.numeric(kcv_43_opt$opt_para[3]),
                 "kcv_63_lossM"=kcv_63_lossM, "kcv_63_loss_truncated"=kcv_63_loss_truncated,"lambda_63" = as.numeric(kcv_63_opt$opt_para[3]),
                 "EV_true" = EV_true,"EV_42" = EV_42,"EV_62" = EV_62, "EV_33" = EV_33,"EV_43" = EV_43,"EV_63" = EV_63
  )
  output
}

save(noise0, file="./rank_misspecification_ev/noise0.Rda")

## noise level 0.5*empirical var
set.seed(20261009)

noise05 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","multiway","funData","MFPCA","SmoothHOOI")) %dorng%{
  sim_data <- sim_data2(L = L_tilde, R= R_tilde, mean_G = mean_G, cov_G = cov_G,
                        E = E, p=200, noise_level=0.5, pattern="random", percent=0.2)
  Mmiss <- sim_data$sim_Mmiss
  Msmooth <- sim_data$sim_Msmooth
  
  ## 10-fold cross-validation--true rank
  kcv_true_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=true_rank_grid, lambda_seq=lambda_seq, k=10,
                             L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_true_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_true_opt$opt_para[1:2]), lambda=as.numeric(kcv_true_opt$opt_para[3]),
                         L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_true_tilde <- MakeIdent(L=kcv_true_res$L, G=kcv_true_res$G, R=kcv_true_res$R)
  
  kcv_true_loss <- loss(tnsr = kcv_true_res$est, true_tnsr=Msmooth@data,
                        L = kcv_true_tilde$L_tilde, true_L = L_tilde,
                        R = kcv_true_tilde$R_tilde, true_R = R_tilde)
  
  
  EV_true <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_true_tilde$L_tilde,
    G_tilde = kcv_true_tilde$G_tilde,
    R_tilde = kcv_true_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--42
  kcv_42_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_42, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_42_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_42_opt$opt_para[1:2]), lambda=as.numeric(kcv_42_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_42_tilde <- MakeIdent(L=kcv_42_res$L, G=kcv_42_res$G, R=kcv_42_res$R)
  
  kcv_42_lossM <- loss(tnsr = kcv_42_res$est, true_tnsr=Msmooth@data)
  
  kcv_42_truncated_res <- rank_truncation(Mmiss, mode="L", true_rank_L=3, true_rank_R=NULL,
                                          L_tilde=kcv_42_tilde$L_tilde, G_tilde=kcv_42_tilde$G_tilde, R_tilde=kcv_42_tilde$R_tilde)
  
  kcv_42_loss_truncated <- loss(tnsr = kcv_42_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_42_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_42_tilde$R_tilde, true_R = R_tilde)
  
  EV_42 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_42_tilde$L_tilde,
    G_tilde = kcv_42_tilde$G_tilde,
    R_tilde = kcv_42_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--62
  kcv_62_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_62, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_62_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_62_opt$opt_para[1:2]), lambda=as.numeric(kcv_62_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_62_tilde <- MakeIdent(L=kcv_62_res$L, G=kcv_62_res$G, R=kcv_62_res$R)
  
  kcv_62_lossM <- loss(tnsr = kcv_62_res$est, true_tnsr=Msmooth@data)
  
  kcv_62_truncated_res <- rank_truncation(Mmiss, mode="L", true_rank_L=3, true_rank_R=NULL,
                                          L_tilde=kcv_62_tilde$L_tilde, G_tilde=kcv_62_tilde$G_tilde, R_tilde=kcv_62_tilde$R_tilde)
  
  kcv_62_loss_truncated <- loss(tnsr = kcv_62_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_62_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_62_tilde$R_tilde, true_R = R_tilde)
  
  EV_62 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_62_tilde$L_tilde,
    G_tilde = kcv_62_tilde$G_tilde,
    R_tilde = kcv_62_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--33
  kcv_33_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_33, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_33_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_33_opt$opt_para[1:2]), lambda=as.numeric(kcv_33_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_33_tilde <- MakeIdent(L=kcv_33_res$L, G=kcv_33_res$G, R=kcv_33_res$R)
  
  kcv_33_lossM <- loss(tnsr = kcv_33_res$est, true_tnsr=Msmooth@data)
  
  kcv_33_truncated_res <- rank_truncation(Mmiss, mode="R", true_rank_L=NULL, true_rank_R=2,
                                          L_tilde=kcv_33_tilde$L_tilde, G_tilde=kcv_33_tilde$G_tilde, R_tilde=kcv_33_tilde$R_tilde)
  
  kcv_33_loss_truncated <- loss(tnsr = kcv_33_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_33_tilde$L_tilde, true_L = L_tilde,
                                R = kcv_33_truncated_res$R_truncated, true_R = R_tilde)
  
  EV_33 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_33_tilde$L_tilde,
    G_tilde = kcv_33_tilde$G_tilde,
    R_tilde = kcv_33_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--43
  kcv_43_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_43, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_43_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_43_opt$opt_para[1:2]), lambda=as.numeric(kcv_43_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_43_tilde <- MakeIdent(L=kcv_43_res$L, G=kcv_43_res$G, R=kcv_43_res$R)
  
  kcv_43_lossM <- loss(tnsr = kcv_43_res$est, true_tnsr=Msmooth@data)
  
  kcv_43_truncated_res <- rank_truncation(Mmiss, mode="both", true_rank_L=3, true_rank_R=2,
                                          L_tilde=kcv_43_tilde$L_tilde, G_tilde=kcv_43_tilde$G_tilde, R_tilde=kcv_43_tilde$R_tilde)
  
  kcv_43_loss_truncated <- loss(tnsr = kcv_43_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_43_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_43_truncated_res$R_truncated, true_R = R_tilde)
  
  EV_43 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_43_tilde$L_tilde,
    G_tilde = kcv_43_tilde$G_tilde,
    R_tilde = kcv_43_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--63
  kcv_63_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_63, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_63_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_63_opt$opt_para[1:2]), lambda=as.numeric(kcv_63_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_63_tilde <- MakeIdent(L=kcv_63_res$L, G=kcv_63_res$G, R=kcv_63_res$R)
  
  kcv_63_lossM <- loss(tnsr = kcv_63_res$est, true_tnsr=Msmooth@data)
  
  kcv_63_truncated_res <- rank_truncation(Mmiss, mode="both", true_rank_L=3, true_rank_R=2,
                                          L_tilde=kcv_63_tilde$L_tilde, G_tilde=kcv_63_tilde$G_tilde, R_tilde=kcv_63_tilde$R_tilde)
  
  kcv_63_loss_truncated <- loss(tnsr = kcv_63_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_63_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_63_truncated_res$R_truncated, true_R = R_tilde)
  
  
  EV_63 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_63_tilde$L_tilde,
    G_tilde = kcv_63_tilde$G_tilde,
    R_tilde = kcv_63_tilde$R_tilde
  )
  
  output <- list("kcv_true_loss"=kcv_true_loss, lambda_true = as.numeric(kcv_true_opt$opt_para[3]),
                 "kcv_42_lossM"=kcv_42_lossM, "kcv_42_loss_truncated"=kcv_42_loss_truncated,"lambda_42" = as.numeric(kcv_42_opt$opt_para[3]),
                 "kcv_62_lossM"=kcv_62_lossM, "kcv_62_loss_truncated"=kcv_62_loss_truncated,"lambda_62" = as.numeric(kcv_62_opt$opt_para[3]),
                 "kcv_33_lossM"=kcv_33_lossM, "kcv_33_loss_truncated"=kcv_33_loss_truncated,"lambda_33" = as.numeric(kcv_33_opt$opt_para[3]),
                 "kcv_43_lossM"=kcv_43_lossM, "kcv_43_loss_truncated"=kcv_43_loss_truncated,"lambda_43" = as.numeric(kcv_43_opt$opt_para[3]),
                 "kcv_63_lossM"=kcv_63_lossM, "kcv_63_loss_truncated"=kcv_63_loss_truncated,"lambda_63" = as.numeric(kcv_63_opt$opt_para[3]),
                 "EV_true" = EV_true,"EV_42" = EV_42,"EV_62" = EV_62, "EV_33" = EV_33,"EV_43" = EV_43,"EV_63" = EV_63
  )
  output
}

save(noise05, file="./rank_misspecification_ev/noise05.Rda")

## noise level 1.5*empirical var
set.seed(20261010)

noise15 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","multiway","funData","MFPCA","SmoothHOOI")) %dorng%{
  sim_data <- sim_data2(L = L_tilde, R= R_tilde, mean_G = mean_G, cov_G = cov_G,
                        E = E, p=200, noise_level=1.5, pattern="random", percent=0.2)
  Mmiss <- sim_data$sim_Mmiss
  Msmooth <- sim_data$sim_Msmooth
  
  ## 10-fold cross-validation--true rank
  kcv_true_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=true_rank_grid, lambda_seq=lambda_seq, k=10,
                             L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_true_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_true_opt$opt_para[1:2]), lambda=as.numeric(kcv_true_opt$opt_para[3]),
                         L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_true_tilde <- MakeIdent(L=kcv_true_res$L, G=kcv_true_res$G, R=kcv_true_res$R)
  
  kcv_true_loss <- loss(tnsr = kcv_true_res$est, true_tnsr=Msmooth@data,
                        L = kcv_true_tilde$L_tilde, true_L = L_tilde,
                        R = kcv_true_tilde$R_tilde, true_R = R_tilde)
  
  
  EV_true <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_true_tilde$L_tilde,
    G_tilde = kcv_true_tilde$G_tilde,
    R_tilde = kcv_true_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--42
  kcv_42_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_42, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_42_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_42_opt$opt_para[1:2]), lambda=as.numeric(kcv_42_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_42_tilde <- MakeIdent(L=kcv_42_res$L, G=kcv_42_res$G, R=kcv_42_res$R)
  
  kcv_42_lossM <- loss(tnsr = kcv_42_res$est, true_tnsr=Msmooth@data)
  
  kcv_42_truncated_res <- rank_truncation(Mmiss, mode="L", true_rank_L=3, true_rank_R=NULL,
                                          L_tilde=kcv_42_tilde$L_tilde, G_tilde=kcv_42_tilde$G_tilde, R_tilde=kcv_42_tilde$R_tilde)
  
  kcv_42_loss_truncated <- loss(tnsr = kcv_42_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_42_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_42_tilde$R_tilde, true_R = R_tilde)
  
  EV_42 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_42_tilde$L_tilde,
    G_tilde = kcv_42_tilde$G_tilde,
    R_tilde = kcv_42_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--62
  kcv_62_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_62, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_62_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_62_opt$opt_para[1:2]), lambda=as.numeric(kcv_62_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_62_tilde <- MakeIdent(L=kcv_62_res$L, G=kcv_62_res$G, R=kcv_62_res$R)
  
  kcv_62_lossM <- loss(tnsr = kcv_62_res$est, true_tnsr=Msmooth@data)
  
  kcv_62_truncated_res <- rank_truncation(Mmiss, mode="L", true_rank_L=3, true_rank_R=NULL,
                                          L_tilde=kcv_62_tilde$L_tilde, G_tilde=kcv_62_tilde$G_tilde, R_tilde=kcv_62_tilde$R_tilde)
  
  kcv_62_loss_truncated <- loss(tnsr = kcv_62_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_62_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_62_tilde$R_tilde, true_R = R_tilde)
  
  EV_62 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_62_tilde$L_tilde,
    G_tilde = kcv_62_tilde$G_tilde,
    R_tilde = kcv_62_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--33
  kcv_33_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_33, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_33_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_33_opt$opt_para[1:2]), lambda=as.numeric(kcv_33_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_33_tilde <- MakeIdent(L=kcv_33_res$L, G=kcv_33_res$G, R=kcv_33_res$R)
  
  kcv_33_lossM <- loss(tnsr = kcv_33_res$est, true_tnsr=Msmooth@data)
  
  kcv_33_truncated_res <- rank_truncation(Mmiss, mode="R", true_rank_L=NULL, true_rank_R=2,
                                          L_tilde=kcv_33_tilde$L_tilde, G_tilde=kcv_33_tilde$G_tilde, R_tilde=kcv_33_tilde$R_tilde)
  
  kcv_33_loss_truncated <- loss(tnsr = kcv_33_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_33_tilde$L_tilde, true_L = L_tilde,
                                R = kcv_33_truncated_res$R_truncated, true_R = R_tilde)
  
  EV_33 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_33_tilde$L_tilde,
    G_tilde = kcv_33_tilde$G_tilde,
    R_tilde = kcv_33_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--43
  kcv_43_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_43, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_43_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_43_opt$opt_para[1:2]), lambda=as.numeric(kcv_43_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_43_tilde <- MakeIdent(L=kcv_43_res$L, G=kcv_43_res$G, R=kcv_43_res$R)
  
  kcv_43_lossM <- loss(tnsr = kcv_43_res$est, true_tnsr=Msmooth@data)
  
  kcv_43_truncated_res <- rank_truncation(Mmiss, mode="both", true_rank_L=3, true_rank_R=2,
                                          L_tilde=kcv_43_tilde$L_tilde, G_tilde=kcv_43_tilde$G_tilde, R_tilde=kcv_43_tilde$R_tilde)
  
  kcv_43_loss_truncated <- loss(tnsr = kcv_43_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_43_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_43_truncated_res$R_truncated, true_R = R_tilde)
  
  EV_43 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_43_tilde$L_tilde,
    G_tilde = kcv_43_tilde$G_tilde,
    R_tilde = kcv_43_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--63
  kcv_63_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_63, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_63_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_63_opt$opt_para[1:2]), lambda=as.numeric(kcv_63_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_63_tilde <- MakeIdent(L=kcv_63_res$L, G=kcv_63_res$G, R=kcv_63_res$R)
  
  kcv_63_lossM <- loss(tnsr = kcv_63_res$est, true_tnsr=Msmooth@data)
  
  kcv_63_truncated_res <- rank_truncation(Mmiss, mode="both", true_rank_L=3, true_rank_R=2,
                                          L_tilde=kcv_63_tilde$L_tilde, G_tilde=kcv_63_tilde$G_tilde, R_tilde=kcv_63_tilde$R_tilde)
  
  kcv_63_loss_truncated <- loss(tnsr = kcv_63_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_63_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_63_truncated_res$R_truncated, true_R = R_tilde)
  
  
  EV_63 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_63_tilde$L_tilde,
    G_tilde = kcv_63_tilde$G_tilde,
    R_tilde = kcv_63_tilde$R_tilde
  )
  
  output <- list("kcv_true_loss"=kcv_true_loss, lambda_true = as.numeric(kcv_true_opt$opt_para[3]),
                 "kcv_42_lossM"=kcv_42_lossM, "kcv_42_loss_truncated"=kcv_42_loss_truncated,"lambda_42" = as.numeric(kcv_42_opt$opt_para[3]),
                 "kcv_62_lossM"=kcv_62_lossM, "kcv_62_loss_truncated"=kcv_62_loss_truncated,"lambda_62" = as.numeric(kcv_62_opt$opt_para[3]),
                 "kcv_33_lossM"=kcv_33_lossM, "kcv_33_loss_truncated"=kcv_33_loss_truncated,"lambda_33" = as.numeric(kcv_33_opt$opt_para[3]),
                 "kcv_43_lossM"=kcv_43_lossM, "kcv_43_loss_truncated"=kcv_43_loss_truncated,"lambda_43" = as.numeric(kcv_43_opt$opt_para[3]),
                 "kcv_63_lossM"=kcv_63_lossM, "kcv_63_loss_truncated"=kcv_63_loss_truncated,"lambda_63" = as.numeric(kcv_63_opt$opt_para[3]),
                 "EV_true" = EV_true,"EV_42" = EV_42,"EV_62" = EV_62, "EV_33" = EV_33,"EV_43" = EV_43,"EV_63" = EV_63
  )
  output
}

save(noise15, file="./rank_misspecification_ev/noise15.Rda")

## noise level 2*empirical var
set.seed(20261011)

noise2 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","multiway","funData","MFPCA","SmoothHOOI")) %dorng%{
  sim_data <- sim_data2(L = L_tilde, R= R_tilde, mean_G = mean_G, cov_G = cov_G,
                        E = E, p=200, noise_level=2, pattern="random", percent=0.2)
  Mmiss <- sim_data$sim_Mmiss
  Msmooth <- sim_data$sim_Msmooth
  
  ## 10-fold cross-validation--true rank
  kcv_true_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=true_rank_grid, lambda_seq=lambda_seq, k=10,
                             L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_true_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_true_opt$opt_para[1:2]), lambda=as.numeric(kcv_true_opt$opt_para[3]),
                         L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_true_tilde <- MakeIdent(L=kcv_true_res$L, G=kcv_true_res$G, R=kcv_true_res$R)
  
  kcv_true_loss <- loss(tnsr = kcv_true_res$est, true_tnsr=Msmooth@data,
                        L = kcv_true_tilde$L_tilde, true_L = L_tilde,
                        R = kcv_true_tilde$R_tilde, true_R = R_tilde)
  
  
  EV_true <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_true_tilde$L_tilde,
    G_tilde = kcv_true_tilde$G_tilde,
    R_tilde = kcv_true_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--42
  kcv_42_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_42, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_42_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_42_opt$opt_para[1:2]), lambda=as.numeric(kcv_42_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_42_tilde <- MakeIdent(L=kcv_42_res$L, G=kcv_42_res$G, R=kcv_42_res$R)
  
  kcv_42_lossM <- loss(tnsr = kcv_42_res$est, true_tnsr=Msmooth@data)
  
  kcv_42_truncated_res <- rank_truncation(Mmiss, mode="L", true_rank_L=3, true_rank_R=NULL,
                                          L_tilde=kcv_42_tilde$L_tilde, G_tilde=kcv_42_tilde$G_tilde, R_tilde=kcv_42_tilde$R_tilde)
  
  kcv_42_loss_truncated <- loss(tnsr = kcv_42_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_42_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_42_tilde$R_tilde, true_R = R_tilde)
  
  EV_42 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_42_tilde$L_tilde,
    G_tilde = kcv_42_tilde$G_tilde,
    R_tilde = kcv_42_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--62
  kcv_62_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_62, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_62_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_62_opt$opt_para[1:2]), lambda=as.numeric(kcv_62_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_62_tilde <- MakeIdent(L=kcv_62_res$L, G=kcv_62_res$G, R=kcv_62_res$R)
  
  kcv_62_lossM <- loss(tnsr = kcv_62_res$est, true_tnsr=Msmooth@data)
  
  kcv_62_truncated_res <- rank_truncation(Mmiss, mode="L", true_rank_L=3, true_rank_R=NULL,
                                          L_tilde=kcv_62_tilde$L_tilde, G_tilde=kcv_62_tilde$G_tilde, R_tilde=kcv_62_tilde$R_tilde)
  
  kcv_62_loss_truncated <- loss(tnsr = kcv_62_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_62_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_62_tilde$R_tilde, true_R = R_tilde)
  
  EV_62 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_62_tilde$L_tilde,
    G_tilde = kcv_62_tilde$G_tilde,
    R_tilde = kcv_62_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--33
  kcv_33_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_33, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_33_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_33_opt$opt_para[1:2]), lambda=as.numeric(kcv_33_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_33_tilde <- MakeIdent(L=kcv_33_res$L, G=kcv_33_res$G, R=kcv_33_res$R)
  
  kcv_33_lossM <- loss(tnsr = kcv_33_res$est, true_tnsr=Msmooth@data)
  
  kcv_33_truncated_res <- rank_truncation(Mmiss, mode="R", true_rank_L=NULL, true_rank_R=2,
                                          L_tilde=kcv_33_tilde$L_tilde, G_tilde=kcv_33_tilde$G_tilde, R_tilde=kcv_33_tilde$R_tilde)
  
  kcv_33_loss_truncated <- loss(tnsr = kcv_33_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_33_tilde$L_tilde, true_L = L_tilde,
                                R = kcv_33_truncated_res$R_truncated, true_R = R_tilde)
  
  EV_33 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_33_tilde$L_tilde,
    G_tilde = kcv_33_tilde$G_tilde,
    R_tilde = kcv_33_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--43
  kcv_43_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_43, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_43_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_43_opt$opt_para[1:2]), lambda=as.numeric(kcv_43_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_43_tilde <- MakeIdent(L=kcv_43_res$L, G=kcv_43_res$G, R=kcv_43_res$R)
  
  kcv_43_lossM <- loss(tnsr = kcv_43_res$est, true_tnsr=Msmooth@data)
  
  kcv_43_truncated_res <- rank_truncation(Mmiss, mode="both", true_rank_L=3, true_rank_R=2,
                                          L_tilde=kcv_43_tilde$L_tilde, G_tilde=kcv_43_tilde$G_tilde, R_tilde=kcv_43_tilde$R_tilde)
  
  kcv_43_loss_truncated <- loss(tnsr = kcv_43_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_43_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_43_truncated_res$R_truncated, true_R = R_tilde)
  
  EV_43 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_43_tilde$L_tilde,
    G_tilde = kcv_43_tilde$G_tilde,
    R_tilde = kcv_43_tilde$R_tilde
  )
  
  
  ## 10-fold cross-validation--63
  kcv_63_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid_63, lambda_seq=lambda_seq, k=10,
                           L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_63_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_63_opt$opt_para[1:2]), lambda=as.numeric(kcv_63_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_63_tilde <- MakeIdent(L=kcv_63_res$L, G=kcv_63_res$G, R=kcv_63_res$R)
  
  kcv_63_lossM <- loss(tnsr = kcv_63_res$est, true_tnsr=Msmooth@data)
  
  kcv_63_truncated_res <- rank_truncation(Mmiss, mode="both", true_rank_L=3, true_rank_R=2,
                                          L_tilde=kcv_63_tilde$L_tilde, G_tilde=kcv_63_tilde$G_tilde, R_tilde=kcv_63_tilde$R_tilde)
  
  kcv_63_loss_truncated <- loss(tnsr = kcv_63_truncated_res$est_truncated, true_tnsr=Msmooth@data,
                                L = kcv_63_truncated_res$L_truncated, true_L = L_tilde,
                                R = kcv_63_truncated_res$R_truncated, true_R = R_tilde)
  
  
  EV_63 <- calc_EV(
    M_denom = Mmiss@data,
    L_tilde = kcv_63_tilde$L_tilde,
    G_tilde = kcv_63_tilde$G_tilde,
    R_tilde = kcv_63_tilde$R_tilde
  )
  
  output <- list("kcv_true_loss"=kcv_true_loss, lambda_true = as.numeric(kcv_true_opt$opt_para[3]),
                 "kcv_42_lossM"=kcv_42_lossM, "kcv_42_loss_truncated"=kcv_42_loss_truncated,"lambda_42" = as.numeric(kcv_42_opt$opt_para[3]),
                 "kcv_62_lossM"=kcv_62_lossM, "kcv_62_loss_truncated"=kcv_62_loss_truncated,"lambda_62" = as.numeric(kcv_62_opt$opt_para[3]),
                 "kcv_33_lossM"=kcv_33_lossM, "kcv_33_loss_truncated"=kcv_33_loss_truncated,"lambda_33" = as.numeric(kcv_33_opt$opt_para[3]),
                 "kcv_43_lossM"=kcv_43_lossM, "kcv_43_loss_truncated"=kcv_43_loss_truncated,"lambda_43" = as.numeric(kcv_43_opt$opt_para[3]),
                 "kcv_63_lossM"=kcv_63_lossM, "kcv_63_loss_truncated"=kcv_63_loss_truncated,"lambda_63" = as.numeric(kcv_63_opt$opt_para[3]),
                 "EV_true" = EV_true,"EV_42" = EV_42,"EV_62" = EV_62, "EV_33" = EV_33,"EV_43" = EV_43,"EV_63" = EV_63
  )
  output
}

save(noise2, file="./rank_misspecification_ev/noise2.Rda")


stopCluster(cl)
