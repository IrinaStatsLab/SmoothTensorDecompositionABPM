# Run FPCA and return the losses 
# tnsr: simulated noisy, incomplete data
# smooth_tnsr: underlying smooth, complete data 
# true L: ground truth of L
# npc, pve, center: arguments required in refund::fpca.sc()
fpca_run <- function(tnsr, npc=NULL, pve=0.99, center=TRUE){
  nmiss_idx <- which(!is.na(tnsr@data))
  b <- tnsr@modes[2]
  p <- tnsr@modes[3]
  
  h <- list()
  outfpca <- list()
  for (i in 1:b){
    h[[i]] <- rTensor::unfold(tnsr[,i,], row_idx = 2, col_idx = 1)@data
    outfpca[[i]] <- refund::fpca.sc(Y = as.matrix(h[[i]]), pve=pve, npc = npc, center=center)
  }
  
  outfpcaYhat <- array(NA, dim=tnsr@modes)
  
  for (i in 1:p){
    for (j in 1:b){
      outfpcaYhat[,j,i] <- outfpca[[j]]$Yhat[i,] # yhat is the estimated, smoothed curve, with missing data imputed
    }
  }

  # take the minimum number of columns among the estimated L matrices
  min_col <- min(sapply(outfpca[1:b], function(x) ncol(x$efunctions)))
  
  # Take average of the estimated L matrices and ensure orthogonality using qr decomposition
  Lfpca <- list()
  for (i in 1:b){
    Lfpca_i <- outfpca[[i]]$efunctions[,1:min_col]
    Lfpca[[i]] <- qr.Q(qr(Lfpca_i))
  }
  
  Lfpca_avg <- Reduce("+", Lfpca) / length(Lfpca)
  Lfpca_avg <- qr.Q(qr(Lfpca_avg))
  
  out=list(est=outfpcaYhat, Lfpca=Lfpca_avg)
}


