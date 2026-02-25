cp_run <- function(tnsr_array, R_seq, const=c("ortsmo", "uncons", "uncons")){
  # ensure all the missing elements are NA, not nan
  nan_idx <- is.nan(tnsr_array)
  tnsr_array[nan_idx] <- NA
  
  nmiss_idx <- which(!is.na(tnsr_array))
  
  if (length(R_seq)==1){
    start <- Sys.time()
    res <- multiway::parafac(X=tnsr_array, nfac=R_seq[1], const=const, output="best", verbose=FALSE)
    end <- Sys.time()
    exe_time <- end - start
    
    res_rescaled <- multiway::rescale(res, mode = "A", newscale = 1/sqrt(nrow(res$A)), absorb = "C")
    
    est <- fitted(res_rescaled)
    Lcp <- res_rescaled$A
    
    ## order L
    norm_vec <- c()
    nfac <- R_seq[1]
    for (r in 1:nfac){
      ar <- res_rescaled$A[,r]
      br <- res_rescaled$B[,r]
      cr <- res_rescaled$C[,r]
      
      norm_vec[r] <- sum(ar^2) * sum(br^2) * sum(cr^2)
    }
    
    Lcp_ordered <- Lcp[, order(norm_vec, decreasing = TRUE)]
    
    
    out=list(est=est, Lcp=Lcp, Lcp_ordered=Lcp_ordered, time = exe_time)
    
  } else if (length(R_seq)>1){
    exp_var = c()
    est = list()
    Lcp = list()
    res_rescaled_list = list()
    for (i in 1:length(R_seq)){
      R_i <- R_seq[i]
      res_i <- multiway::parafac(X=tnsr_array, nfac=R_i, const=const, output="best",verbose = FALSE)
      
      res_rescaled_i <- multiway::rescale(res_i, mode = "A", newscale = 1/sqrt(nrow(res_i$A)), absorb = "C")
      
      est_i <- fitted(res_rescaled_i)
      Lcp_i <- res_rescaled_i$A
      
      res_rescaled_list[[i]] = res_rescaled_i
      exp_var[i] <- sum((est_i[nmiss_idx])^2)/sum((tnsr_array[nmiss_idx])^2)
      est[[i]] <- est_i
      Lcp[[i]] <- Lcp_i
    }
    best_R_idx <- which.max(exp_var)
    best_R <- R_seq[best_R_idx]
    best_est <- est[[best_R_idx]]
    best_Lcp <- Lcp[[best_R_idx]]
    best_res_rescaled <- res_rescaled_list[[best_R_idx]]
    
    ## order L
    norm_vec <- c()
    nfac <- best_R
    for (r in 1:nfac){
      ar <- best_res_rescaled$A[,r]
      br <- best_res_rescaled$B[,r]
      cr <- best_res_rescaled$C[,r]
      
      norm_vec[r] <- sum(ar^2) * sum(br^2) * sum(cr^2)
    }
    
    best_Lcp_ordered <- best_Lcp[, order(norm_vec, decreasing = TRUE)]
    
    out=list(best_R = best_R, est=best_est, Lcp=best_Lcp, Lcp_ordered = best_Lcp_ordered)
    
  }
}


