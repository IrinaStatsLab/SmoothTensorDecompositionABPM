library(Rcpp)
library(RcppArmadillo)
library(doParallel)
library(doRNG)
library(rTensor)
library(MASS)
library(refund)
library(SmoothHOOI)

## Should run this on HPC 

nworkers <- detectCores()
cl <- makePSOCKcluster(nworkers)
registerDoParallel(cl)

# Results generated from running SmoothHOOI algorithm on real data, with hyperparameter optimization and identifiability correction
load("/home/lyqian/ABPM_new/AlgorithmResIdent.Rda") 

N <- 100
p <- 200
noise_level <- 1
D2 <- SecDiffMat(24)
para_grid <- as.matrix(expand.grid(r1<-c(3), r2<-c(2), lambda <- seq(1,50,by=1)))
rank_grid <- as.matrix(expand.grid(r1<-c(3), r2<-c(2)))
lambda_seq <- seq(1,50,by=1)


## no missingness
set.seed(98765)

fixr_miss0 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","SmoothHOOI")) %dorng%{
  ## Generate a simulated data
  sim_data <- simdata_generator(L_tilde, G_tilde, R_tilde, E, p=200, noise_level=1, pattern="random", percent=0)
  Mmiss <- sim_data$sim_Mmiss
  Msmooth <- sim_data$sim_Msmooth
  
  ## Oracle
  oracle_opt <- oracle(tnsr=Mmiss@data, smooth_tnsr=Msmooth@data, rank_grid=rank_grid, lambda_seq=lambda_seq,
                       init=0, D = D2, max_iter = 500, tol = 0.1, L0 = NULL)
  
  oracle_res <- mglram(Mmiss@data, ranks=as.numeric(oracle_opt$opt_para[1:2]), lambda=as.numeric(oracle_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  oracle_tilde <- MakeIdent(L=oracle_res$L, G=oracle_res$G, R=oracle_res$R)
  
  oracle_loss <- loss(tnsr = oracle_res$est, smooth_tnsr=Msmooth@data,
                      L = oracle_tilde$L_tilde, true_L = L_tilde,
                      R = oracle_tilde$R_tilde, true_R = R_tilde)
  
  ## 10-fold cross-validation
  kcv_opt <- kcv(tnsr=Mmiss@data, rank_grid=rank_grid, lambda_seq=lambda_seq, k=10,
                 L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_opt$opt_para[1:2]), lambda=as.numeric(kcv_opt$opt_para[3]),
                    L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  ## Identifiability Correction
  kcv_tilde <- MakeIdent(L=kcv_res$L, G=kcv_res$G, R=kcv_res$R)
  
  ## Calculate losses 
  kcv_loss <- loss(tnsr = kcv_res$est, smooth_tnsr=Msmooth@data,
                   L = kcv_tilde$L_tilde, true_L = L_tilde,
                   R = kcv_tilde$R_tilde, true_R = R_tilde)
  
  ## FPCA, fix number of principal components to 3
  fpca_loss <- fpca_res(Mmiss, smooth_tnsr=Msmooth, true_L=L_tilde, npc=3, center=FALSE)
  
  output <- list("oracle_loss"=oracle_loss, "kcv_loss"=kcv_loss, "fpca_loss"=fpca_loss)
  output
}

# save(fixr_miss0, file="/home/lyqian/ABPM_new/cf3_res_miss0.Rda")

## missing rate 10% 
set.seed(97531)

fixr_miss10 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","SmoothHOOI")) %dorng%{
  sim_data <- simdata_generator(L_tilde, G_tilde, R_tilde, E, p=200, noise_level=1, pattern="random", percent=0.1)
  Mmiss <- sim_data$sim_Mmiss
  Msmooth <- sim_data$sim_Msmooth
  
  oracle_opt <- oracle(tnsr=Mmiss@data, smooth_tnsr=Msmooth@data, rank_grid=rank_grid, lambda_seq=lambda_seq,
                       init=0, D = D2, max_iter = 500, tol = 0.1, L0 = NULL)
  
  oracle_res <- mglram(Mmiss@data, ranks=as.numeric(oracle_opt$opt_para[1:2]), lambda=as.numeric(oracle_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  oracle_tilde <- MakeIdent(L=oracle_res$L, G=oracle_res$G, R=oracle_res$R)
  
  oracle_loss <- loss(tnsr = oracle_res$est, smooth_tnsr=Msmooth@data,
                      L = oracle_tilde$L_tilde, true_L = L_tilde,
                      R = oracle_tilde$R_tilde, true_R = R_tilde)
  
  kcv_opt <- kcv(tnsr=Mmiss@data, rank_grid=rank_grid, lambda_seq=lambda_seq, k=10,
                 L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_opt$opt_para[1:2]), lambda=as.numeric(kcv_opt$opt_para[3]),
                    L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_tilde <- MakeIdent(L=kcv_res$L, G=kcv_res$G, R=kcv_res$R)
  
  kcv_loss <- loss(tnsr = kcv_res$est, smooth_tnsr=Msmooth@data,
                   L = kcv_tilde$L_tilde, true_L = L_tilde,
                   R = kcv_tilde$R_tilde, true_R = R_tilde)
  
  fpca_loss <- fpca_res(Mmiss, smooth_tnsr=Msmooth, true_L=L_tilde, npc=3, center=FALSE)
  
  output <- list("oracle_loss"=oracle_loss, "kcv_loss"=kcv_loss, "fpca_loss"=fpca_loss)
  output
}

# save(fixr_miss10, file="/home/lyqian/ABPM_new/cf3_res_miss10.Rda")

## missing rate 20% 
set.seed(86420)

fixr_miss20 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","SmoothHOOI")) %dorng%{
  sim_data <- simdata_generator(L_tilde, G_tilde, R_tilde, E, p=200, noise_level=1, pattern="random", percent=0.2)
  Mmiss <- sim_data$sim_Mmiss
  Msmooth <- sim_data$sim_Msmooth
  
  oracle_opt <- oracle(tnsr=Mmiss@data, smooth_tnsr=Msmooth@data, rank_grid=rank_grid, lambda_seq=lambda_seq,
                       init=0, D = D2, max_iter = 500, tol = 0.1, L0 = NULL)
  
  oracle_res <- mglram(Mmiss@data, ranks=as.numeric(oracle_opt$opt_para[1:2]), lambda=as.numeric(oracle_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  oracle_tilde <- MakeIdent(L=oracle_res$L, G=oracle_res$G, R=oracle_res$R)
  
  oracle_loss <- loss(tnsr = oracle_res$est, smooth_tnsr=Msmooth@data,
                      L = oracle_tilde$L_tilde, true_L = L_tilde,
                      R = oracle_tilde$R_tilde, true_R = R_tilde)
  
  kcv_opt <- kcv(tnsr=Mmiss@data, rank_grid=rank_grid, lambda_seq=lambda_seq, k=10,
                 L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_opt$opt_para[1:2]), lambda=as.numeric(kcv_opt$opt_para[3]),
                    L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_tilde <- MakeIdent(L=kcv_res$L, G=kcv_res$G, R=kcv_res$R)
  
  kcv_loss <- loss(tnsr = kcv_res$est, smooth_tnsr=Msmooth@data,
                   L = kcv_tilde$L_tilde, true_L = L_tilde,
                   R = kcv_tilde$R_tilde, true_R = R_tilde)
  
  fpca_loss <- fpca_res(Mmiss, smooth_tnsr=Msmooth, true_L=L_tilde, npc=3, center=FALSE)
  
  output <- list("oracle_loss"=oracle_loss, "kcv_loss"=kcv_loss, "fpca_loss"=fpca_loss)
  output
}

# save(fixr_miss20, file="/home/lyqian/ABPM_new/cf3_res_miss20.Rda")

## missing rate 50% 
set.seed(15151)

fixr_miss50 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","SmoothHOOI")) %dorng%{
  sim_data <- simdata_generator(L_tilde, G_tilde, R_tilde, E, p=200, noise_level=1, pattern="random", percent=0.5)
  Mmiss <- sim_data$sim_Mmiss
  Msmooth <- sim_data$sim_Msmooth
  
  oracle_opt <- oracle(tnsr=Mmiss@data, smooth_tnsr=Msmooth@data, rank_grid=rank_grid, lambda_seq=lambda_seq,
                       init=0, D = D2, max_iter = 500, tol = 0.1, L0 = NULL)
  
  oracle_res <- mglram(Mmiss@data, ranks=as.numeric(oracle_opt$opt_para[1:2]), lambda=as.numeric(oracle_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  oracle_tilde <- MakeIdent(L=oracle_res$L, G=oracle_res$G, R=oracle_res$R)
  
  oracle_loss <- loss(tnsr = oracle_res$est, smooth_tnsr=Msmooth@data,
                      L = oracle_tilde$L_tilde, true_L = L_tilde,
                      R = oracle_tilde$R_tilde, true_R = R_tilde)
  
  kcv_opt <- kcv(tnsr=Mmiss@data, rank_grid=rank_grid, lambda_seq=lambda_seq, k=10,
                 L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_opt$opt_para[1:2]), lambda=as.numeric(kcv_opt$opt_para[3]),
                    L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_tilde <- MakeIdent(L=kcv_res$L, G=kcv_res$G, R=kcv_res$R)
  
  kcv_loss <- loss(tnsr = kcv_res$est, smooth_tnsr=Msmooth@data,
                   L = kcv_tilde$L_tilde, true_L = L_tilde,
                   R = kcv_tilde$R_tilde, true_R = R_tilde)
  
  fpca_loss <- fpca_res(Mmiss, smooth_tnsr=Msmooth, true_L=L_tilde, npc=3, center=FALSE)
  
  output <- list("oracle_loss"=oracle_loss, "kcv_loss"=kcv_loss, "fpca_loss"=fpca_loss)
  output
}

# save(fixr_miss50, file="/home/lyqian/ABPM_new/cf3_res_miss50.Rda")

stopCluster(cl)

