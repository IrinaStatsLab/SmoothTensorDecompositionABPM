load("./missing_normalized_3.Rda")

missing_normalized_3@data[is.nan(missing_normalized_3@data)] <- NA

## CP3 with orthogonality and smoothness
plots <- list()
for (i in 1:10){
  res <- multiway::parafac(X=missing_normalized_3@data, nfac=c(3), const=c("ortsmo", "uncons", "uncons"), output="best")
  
  res_rescaled <- multiway::rescale(res, mode = "A", newscale = 1/sqrt(nrow(res$A)), absorb = "C")
  
  dataL = data.frame(loading = c(res_rescaled$A[, 1], res_rescaled$A[, 2], res_rescaled$A[, 3]),
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
  
  plots[[i]] <- pL
}

## Align plots in a 5*2 grid
library(gridExtra)
combined <- grid.arrange(grobs = plots, ncol = 2, nrow = 5)
combined

# ggsave(plot=combined, filename="CP3ortsmo_L_10runs.pdf", dpi=600, width=12, height=10)

## CP3 with smoothness
plots <- list()
for (i in 1:10){
  res <- multiway::parafac(X=missing_normalized_3@data, nfac=c(3), const=c("smooth", "uncons", "uncons"), output="best")
  
  res_rescaled <- multiway::rescale(res, mode = "A", newscale = 1/sqrt(nrow(res$A)), absorb = "C")
  
  Q <- qr.Q(qr(res_rescaled$A))
  
  dataL = data.frame(loading = c(Q[, 1], Q[, 2], Q[, 3]),
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
  
  plots[[i]] <- pL
}

## Align plots in a 5*2 grid
library(gridExtra)
combined <- grid.arrange(grobs = plots, ncol = 2, nrow = 5)
combined

#ggsave(plot=combined, filename="CP3smooth_L_10runs.pdf", dpi=600, width=12, height=10)

## CP6 with orthogonality and smoothness
plots <- list()
for (i in 1:10){
  res <- multiway::parafac(X=missing_normalized_3@data, nfac=c(6), const=c("ortsmo", "uncons", "uncons"), output="best")
  
  res_rescaled <- multiway::rescale(res, mode = "A", newscale = 1/sqrt(nrow(res$A)), absorb = "C")
  
  dataL = data.frame(loading = c(res_rescaled$A[, 1], res_rescaled$A[, 2], res_rescaled$A[, 3], res_rescaled$A[, 4], res_rescaled$A[, 5], res_rescaled$A[,6]),
                     component = c(rep("1st", 24), rep("2nd", 24), rep("3rd", 24), rep("4th", 24), rep("5th", 24), rep("6th", 24)),
                     hour = rep(12:35, 6))
  
  pL = dataL %>%
    ggplot(aes(x = hour, y = loading)) + geom_line(alpha = 0.5) + geom_point() + geom_hline(yintercept = 0, color = "red") +
    scale_x_continuous(breaks = c(12,  18, 24, 30,  36), labels = c(12,  18,  0,  6, 12)) +
    xlab("Hour of the day") + 
    facet_grid(~ component) + 
    ylab("Coefficient") +
    theme(text = element_text(size = 12)) 
  print(pL)
  
  plots[[i]] <- pL
}

## Align plots in a 5*2 grid
library(gridExtra)
combined <- grid.arrange(grobs = plots, ncol = 2, nrow = 5)
combined

#ggsave(plot=combined, filename="CP6ortsmo_L_10runs.pdf", dpi=600, width=15, height=10)

## CP6 with smoothness
plots <- list()
for (i in 1:10){
  res <- multiway::parafac(X=missing_normalized_3@data, nfac=c(6), const=c("smooth", "uncons", "uncons"), output="best")
  
  res_rescaled <- multiway::rescale(res, mode = "A", newscale = 1/sqrt(nrow(res$A)), absorb = "C")
  
  Q <- qr.Q(qr(res_rescaled$A))
  
  dataL = data.frame(loading = c(Q[, 1], Q[, 2], Q[, 3], Q[,4], Q[,5], Q[,6]),
                     component = c(rep("1st", 24), rep("2nd", 24), rep("3rd", 24), rep("4th", 24), rep("5th", 24), rep("6th", 24)),
                     hour = rep(12:35, 6))
  
  pL = dataL %>%
    ggplot(aes(x = hour, y = loading)) + geom_line(alpha = 0.5) + geom_point() + geom_hline(yintercept = 0, color = "red") +
    scale_x_continuous(breaks = c(12, 15, 18, 21, 24, 27, 30, 33, 36), labels = c(12, 15, 18, 21, 0, 3, 6, 9, 12)) +
    xlab("Hour of the day") + 
    facet_grid(~ component) + 
    ylab("Coefficient") +
    theme(text = element_text(size = 12)) 
  print(pL)
  
  plots[[i]] <- pL
}

## Align plots in a 5*2 grid
library(gridExtra)
combined <- grid.arrange(grobs = plots, ncol = 2, nrow = 5)
combined

#ggsave(plot=combined, filename="CP6smooth_L_10runs.pdf", dpi=600, width=15, height=10)

