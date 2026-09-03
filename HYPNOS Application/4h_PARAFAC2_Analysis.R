library(multiway)
library(dplyr)
library(ggplot2)

load("./missing_normalized_3.Rda")

missing_normalized_3@data[is.nan(missing_normalized_3@data)] <- NA


load("./AlgorithmResIdent.Rda")

## nfac=3, ortsmo

set.seed(123)
res <- multiway::parafac2(X=missing_normalized_3@data, nfac=c(3), const=c("ortsmo", "uncons", "uncons"), output="best")

M <- Reduce(`+`, lapply(res$A, function(A) abs(crossprod(A, L_tilde))))

perm <- integer(ncol(L_tilde))
available_col <- seq_len(ncol(M))
for (j in seq_len(ncol(L_tilde))) {
  k <- available_col[which.max(M[available_col, j])]
  perm[j] <- k
  available_col <- setdiff(available_col, k)
}
perm   

align_to_ref <- function(A, Ref){
  A_perm <- A[, perm, drop = FALSE]
  sgn <- sign(diag(crossprod(A_perm, Ref)))
  sgn[sgn == 0] <- 1
  sweep(A_perm, 2, sgn, "*")
}

df_p2 <- do.call(rbind, lapply(seq_along(res$A), function(i){
  Aa <- align_to_ref(res$A[[i]], L_tilde)
  data.frame(hour = rep(12:35, 3), loading = as.vector(Aa),
             component = rep(c("1st", "2nd", "3rd"), each = 24), id = i)
}))

df_L <- data.frame(hour = rep(12:35, 3), loading = as.vector(L_tilde),
                   component = rep(c("1st", "2nd", "3rd"), each = 24))

ggplot(df_p2, aes(hour, loading)) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = .3) +
  geom_line(aes(group = id), colour = "grey70", alpha = .2, linewidth = .3) +
  geom_line(data = df_L, colour = "black", linewidth = 1.1) +
  geom_point(data = df_L, colour = "black", size = 1.3) +
  facet_grid(~ component) +
  scale_x_continuous(breaks = seq(12, 36, 3),
                     labels = c(12, 15, 18, 21, 0, 3, 6, 9, 12)) +
  labs(x = "Hour of the day", y = "Coefficient") +
  theme_bw(base_size = 12)

#ggsave("./L_comp_PARAFAC2_ortsmo.pdf", dpi=600, width=8, height=4)

## nfac=3, smooth
set.seed(111)
res <- multiway::parafac2(X=missing_normalized_3@data, nfac=c(3), const=c("smooth", "uncons", "uncons"), output="best")


for (i in 1:207){
  
  L <- res$A[[i]]
  
  dataL = data.frame(loading = c(L[, 1], L[, 2], L[, 3]),
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
  
}

