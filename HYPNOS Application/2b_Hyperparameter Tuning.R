load("./missing_normalized_3.Rda")

library(SmoothHOOI)
library(ggplot2)

# Second Order Difference matrix 
D2 <- SecDiffMat(24)

# Record index of missingness in the tensor
nmiss_idx <- which(!is.na(missing_normalized_3@data))

# Run k-fold cross validation, will run for about 3 minutes
set.seed(678678)

start_time <- Sys.time()

kcv_res <- kcv(tnsr=missing_normalized_3@data, rank_grid=as.matrix(expand.grid(r1<-seq(3,6,by=1), r2<-c(2,3))), lambda_seq=seq(1,20,by=1),
               k=5, L0=NULL, D=D2, tol=0.01, max_iter=500, init=0)

end_time <- Sys.time()

exec_time <- end_time - start_time
print(exec_time)

kcv_res$opt_para # 6, 3, 12

kcv_res$MSE_mat

kcv_res$SE_mat 

# Run our algorithm with the "optimal" hyperparameters
res_opt <-  mglram(tnsr = missing_normalized_3@data, ranks = c(6, 3), init=0, D = D2,
                   lambda = 12, max_iter = 500, tol = 1e-5, L0 = NULL)
res_opt$conv

res_opt$L
res_opt$R 

# Cumulative explained variability when r1 increases
tilde <- MakeIdent(res_opt$L, res_opt$G, res_opt$R)
L_tilde <- tilde$L_tilde
G_tilde <- tilde$G_tilde
R_tilde <- tilde$R_tilde

cum_var_L <- rep(NA, 6)
for (i in 1:6){
  comp <- array(NA, dim(missing_normalized_3@data))
  for (j in 1:207){
    comp[ , , j] <- matrix(L_tilde[, 1:i], ncol=i) %*% matrix(G_tilde[1:i, , j], nrow=i) %*% t(R_tilde)
  }
  cum_var_L[i] <- 1 - sum((missing_normalized_3@data[nmiss_idx]-comp[nmiss_idx])^2)/sum(missing_normalized_3@data[nmiss_idx]^2)
}

cum_var_L # 0.4508818 0.6144494 0.6723581 0.6970985 0.7201655 0.7276182

# Separate explained variability for different L ranks
cum_var_L0 <- c(0, cum_var_L[1:5])
sep_var_L <- cum_var_L - cum_var_L0

sep_var_L # 0.45088178 0.16356766 0.05790869 0.02474033 0.02306705 0.00745271

# Cumulative explained variability when r2 increases
cum_var_R <- rep(NA, 3)
for (i in 1:3){
  comp <- array(NA, dim(missing_normalized_3@data))
  for (j in 1:207){
    comp[ , , j] <- L_tilde %*% matrix(G_tilde[,1:i, j], ncol=i) %*% matrix(t(R_tilde[, 1:i]), nrow=i)
  }
  #cum_var_R[i] <- sum(comp[nmiss_idx]^2)/sum(missing_normalized_3@data[nmiss_idx]^2)
  cum_var_R[i] <- 1 - sum((comp[nmiss_idx]-missing_normalized_3@data[nmiss_idx])^2)/sum(missing_normalized_3@data[nmiss_idx]^2)
}

cum_var_R # 0.4275989 0.6483961 0.7276182

# Separate explained variability for different R ranks
cum_var_R0 <- c(0, cum_var_R[1:2])
sep_var_R <- cum_var_R - cum_var_R0

sep_var_R # 0.42759894 0.22079714 0.07922215

# Plot the explained variability
var_data_L <- data.frame(
  component = factor(seq_along(sep_var_L), levels = seq_along(sep_var_L)), 
  value = sort(sep_var_L, decreasing = TRUE) 
)

var_data_L$percentage <- paste0(round(var_data_L$value * 100, 2), "%")

ggplot(var_data_L, aes(x = component, y = value)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(aes(label = percentage), vjust = -0.5, size = 5) +
  geom_line(aes(group = 1), color = "black", linewidth = 0.5) +
  geom_point(color = "black", size = 1) +
  #labs(y = "Explained Variability", x="L Component", title="Explained Variability of the Top 6 Components") +
  labs(y = "Explained Variability", x="L Component", title="") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15), 
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) 

# ggsave("var_L_upodated.pdf", dpi=600, width=6, height=4)

# The 4th - 6th components explain little variability, can be removed

var_data_R <- data.frame(
  component = factor(seq_along(sep_var_R), levels = seq_along(sep_var_R)), 
  value = sort(sep_var_R, decreasing = TRUE) 
)

var_data_R$percentage <- paste0(round(var_data_R$value * 100, 2), "%")

ggplot(var_data_R, aes(x = component, y = value)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 0.4) +
  geom_text(aes(label = percentage), vjust = -0.5, size = 5) +
  geom_line(aes(group = 1), color = "black", linewidth = 0.5) +
  geom_point(color = "black", size = 1) +
  #labs(title = "Explained Variability of the 3 Components", x = "R Component", y = "Explained Variability") +
  labs(title = "", x = "R Component", y = "Explained Variability") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15), 
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# ggsave("var_R_updated.pdf", dpi=600, width=6, height=4)

# The 3rd components explain variability can be removed

# The final rank will be r1=3, r2=2, and the corresponding lambda will be 4

res_adj <-  mglram(tnsr = missing_normalized_3@data, ranks = c(3, 2), init=0, D = D2,
                   lambda = 4, max_iter = 500, tol = 1e-5, L0 = NULL)
res_adj$conv

1 - sum((missing_normalized_3@data[nmiss_idx] - res_adj$est[nmiss_idx])^2)/sum(missing_normalized_3@data[nmiss_idx]^2) # Explained variability of the final model
# 0.6059719

