library(funData)
library(MFPCA)
library(ggplot2)

load("./missing_normalized_3.Rda")

matrix_DBP <- t(missing_normalized_3@data[,1,])
matrix_SBP <- t(missing_normalized_3@data[,2,])
matrix_HR  <- t(missing_normalized_3@data[,3,])

time_grid <- seq(1, 24, by = 1)
fd_DBP <- funData(argvals = list(time_grid), X = matrix_DBP)
fd_SBP <- funData(argvals = list(time_grid), X = matrix_SBP)
fd_HR <- funData(argvals = list(time_grid), X = matrix_HR)

mfd_all <- multiFunData(fd_DBP, fd_SBP, fd_HR)
mfd_all

uniExpansions <- list(list(type = "uFPCA"),
                      list(type = "uFPCA"),
                      list(type = "uFPCA"))

MFPCAres3 <- MFPCA(mfd_all, M = 3, uniExpansions = uniExpansions, fit=TRUE)

plot(MFPCAres3)

MFPCAeigenfunc_to_df <- function(MFPCAeigenfunc, var_names){
  
  df_list <- list()
  
  for (i in 1:length(MFPCAeigenfunc)) {
    data_mat <- MFPCAeigenfunc[[i]]@X
    argvals  <- MFPCAeigenfunc[[i]]@argvals[[1]]
    n_obs    <- nrow(data_mat)
    
    temp_df <- data.frame(
      argvals = rep(argvals, n_obs),
      Value = as.vector(t(data_mat)), 
      Index = rep(1:n_obs, each = length(argvals)), 
      Variable = var_names[i]
    )
    
    df_list[[i]] <- temp_df
  }
  final_df <- do.call(rbind, df_list)
  return(final_df)
}

df_for_plot <- MFPCAeigenfunc_to_df(MFPCAres3$functions, var_names=c("DBP","SBP","HR"))

df_for_plot$Variable <- factor(df_for_plot$Variable, levels = c("DBP", "SBP", "HR"))


my_colors <- c("#E69F00", "#009E73", "#0072B2")
ggplot(df_for_plot, aes(x = argvals, y = Value, color = as.factor(Index))) +
  geom_hline(yintercept = 0, color = "gray", size = 0.8) +
  geom_line(linewidth = 1.2) +
  scale_x_continuous(breaks = c(1, 4, 7, 10, 13, 16, 19, 22, 25), labels = c(12, 15, 18, 21, 0, 3, 6, 9, 12)) +
  facet_wrap(~ Variable, scales = "fixed", ncol = 3) +
  scale_color_manual(values = my_colors) +
  theme_bw() +
  theme(
    legend.position = "bottom",       
    strip.text = element_text(size = 12), 
    axis.title = element_text(size = 11)
  ) +
  labs(title = "Multivariate Functional Principal Components",
       y = "Eigenfunction Value",
       color = "Principal Components (PCs)") +
  guides(color = guide_legend(nrow = 1))

#ggsave(filename="MFPCA_PCs.pdf", dpi=600, width=8, height=3)

