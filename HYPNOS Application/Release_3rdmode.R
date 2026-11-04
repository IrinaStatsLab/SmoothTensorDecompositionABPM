library(rTensor)
library(SmoothHOOI)
library(dplyr)
library(ggplot2)

load("./missing_normalized_3.Rda")

D2 <- SecDiffMat(24)
res <- mglram_freeC(tnsr = missing_normalized_3@data, ranks = c(3, 2, 6), lambda = 4,
                    D = D2, tol = 1e-5, max_iter = 500, init = 0)

res$obj_vec

dataL = data.frame(loading = c(res$L[, 1], res$L[, 2], res$L[, 3]),
                   component = c(rep("1st", 24), rep("2nd", 24), rep("3rd", 24)),
                   hour = rep(12:35, 3))

pL = dataL %>%
  ggplot(aes(x = hour, y = loading)) + geom_line(alpha = 0.5) + geom_point() + geom_hline(yintercept = 0, color = "red") +
  scale_x_continuous(breaks = c(12, 15, 18, 21, 24, 27, 30, 33, 36), labels = c(12, 15, 18, 21, 0, 3, 6, 9, 12)) +
  xlab("Hour of the day") + 
  facet_grid(~ component) + 
  ylab("Coefficient") +
  theme(text = element_text(size = 12)) 
print(pL)

res_SmoothHOOI <- mglram(tnsr = missing_normalized_3@data, ranks = c(3, 2), init=0, D = D2,
                         lambda = 4, max_iter = 500, tol = 1e-5, L0 = NULL)

loss(tnsr=res_SmoothHOOI$est, true_tnsr=res$est, 
     L = res_SmoothHOOI$L, true_L = res$L)

res <- mglram_freeC(tnsr = missing_normalized_3@data, ranks = c(3, 2, 3), lambda = 4,
                    D = D2, tol = 1e-5, max_iter = 500, init = 0)

dataL = data.frame(loading = c(res$L[, 1], res$L[, 2], res$L[, 3]),
                   component = c(rep("1st", 24), rep("2nd", 24), rep("3rd", 24)),
                   hour = rep(12:35, 3))

pL = dataL %>%
  ggplot(aes(x = hour, y = loading)) + geom_line(alpha = 0.5) + geom_point() + geom_hline(yintercept = 0, color = "red") +
  scale_x_continuous(breaks = c(12, 15, 18, 21, 24, 27, 30, 33, 36), labels = c(12, 15, 18, 21, 0, 3, 6, 9, 12)) +
  xlab("Hour of the day") + 
  facet_grid(~ component) + 
  ylab("Coefficient") +
  theme(text = element_text(size = 12)) 
print(pL)

#ggsave("./compress_3rdmode_L.pdf", pL, width = 9, height = 3, dpi = 600)

loss(tnsr=res_SmoothHOOI$est, true_tnsr=res$est, 
     L = res_SmoothHOOI$L, true_L = res$L)

nmiss_idx <- which(!is.na(missing_normalized_3@data))
1 - sum((missing_normalized_3@data[nmiss_idx]-res$est[nmiss_idx])^2)/sum(missing_normalized_3@data[nmiss_idx]^2)



