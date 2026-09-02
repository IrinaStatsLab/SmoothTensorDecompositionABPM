parafac2_run <- function(tnsr_array, R_seq, const=c("ortsmo", "uncons", "uncons")){
  # ensure all the missing elements are NA, not nan
  nan_idx <- is.nan(tnsr_array)
  tnsr_array[nan_idx] <- NA
  
  nmiss_idx <- which(!is.na(tnsr_array))
  
  if (length(R_seq)==1){
    res <- multiway::parafac2(X=tnsr_array, nfac=R_seq[1], const=const, output="best", verbose=FALSE)
    
    est <- fitted(res)
    
    out=list(est=est, Lparafac2_list=res$A)
    return(out)
    
  } else if (length(R_seq)>1){
    exp_var = c()
    est = list()
    for (i in 1:length(R_seq)){
      R_i <- R_seq[i]
      res_i <- multiway::parafac2(X=tnsr_array, nfac=R_i, const=const, output="best",verbose = FALSE)
    
      est_i <- fitted(res_i)
      
      exp_var[i] <- sum((est_i[nmiss_idx])^2)/sum((tnsr_array[nmiss_idx])^2)
      #exp_var[i] <- 1 - sum((tnsr_array[nmiss_idx] - est_i[nmiss_idx])^2) /sum(tnsr_array[nmiss_idx]^2)
      est[[i]] <- est_i
    }
    best_R_idx <- which.max(exp_var)
    best_R <- R_seq[best_R_idx]
    best_est <- est[[best_R_idx]]
    
    out=list(best_R = best_R, est=best_est)
  }
}