library(reshape2)
library(dplyr)
library(ggplot2)
library(SmoothHOOI)

## Read the data from original files 
data_folder <- "./PRSA_Data" 

file_paths <- list.files(path = data_folder, pattern = "*.csv", full.names = TRUE)

file_names <- basename(file_paths)

station_names <- sapply(strsplit(file_names, "_"), function(x) x[3])

## Only keep the air pollutant variables
vars_to_keep <- c("PM2.5", "PM10", "SO2", "NO2", "CO", "O3")

## Take averages across different years 
avg_across <- function(filepath) {
  df <- read.csv(filepath) 
  
  agg_by_day <- df %>%
    dplyr::select(month, day, all_of(vars_to_keep)) %>% 
    group_by(month, day) %>%
    summarise(across(everything(), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>% 
    arrange(month, day) # order the time from Jan 1 to Dec 31
  
  return(as.matrix(agg_by_day[, vars_to_keep]))
}

list_by_station <- lapply(file_paths, avg_across)

air_tensor <- simplify2array(list_by_station) # this is the tensor of dimension (366, 6, 12)

## Standardize the data wrt each variable
dimnames(air_tensor) <- list(
  Time = 1:366,
  Variables = vars_to_keep,
  Stations = station_names
)

air_tensor_std <- array(NA, dim=dim(air_tensor), dimnames = dimnames(air_tensor))

for (var in vars_to_keep) {
  
  var_slice <- air_tensor[, var, ]
  
  slice_mean <- mean(var_slice, na.rm = TRUE)
  slice_sd <- sd(var_slice, na.rm = TRUE)
  
  air_tensor_std[, var, ] <- (var_slice - slice_mean) / slice_sd
  
  cat(sprintf("  %s: Mean=%.2f, SD=%.2f \n", var, slice_mean, slice_sd))
}

# PM2.5: Mean=79.71, SD=35.73 
# PM10: Mean=104.52, SD=40.59 
# SO2: Mean=15.80, SD=11.84 
# NO2: Mean=50.55, SD=18.91 
# CO: Mean=1232.44, SD=633.19 
# O3: Mean=57.18, SD=31.83 

sum(is.na(air_tensor_std)) ## no missing data exist

## Overview of the data 
for (j in 1:12){
  for (i in 1:6){
    plot(air_tensor_std[,i,j])
  }
}


## Make the second difference matrix 366*366
D2 <- SecDiffMat(366)

## Check residual autocorrelation
res_ac <- mglram(tnsr = air_tensor_std, ranks = c(3,3), init=0, D = D2,
              lambda = 100000, max_iter = 500, tol = 1e-5, L0 = NULL)

pdf("PM25_ResidualAutocorrelation.pdf", width = 12, height = 8)
for(i in 1:12){
  if(sum(is.na(air_tensor_std[,1,i])) == 0){
    my_residuals <-air_tensor_std[,1,i] - res_ac$est[,1,i]
    
    acf(my_residuals, 
        na.action = na.pass, main = paste("Monitoring Site ", i, ": PM2.5 Residual Autocorrelation", sep=""),
        lag.max = 120) 
  }
}
dev.off()

pdf("PM10_ResidualAutocorrelation.pdf", width = 12, height = 8)
for(i in 1:12){
  if(sum(is.na(air_tensor_std[,2,i])) == 0){
    my_residuals <-air_tensor_std[,2,i] - res_ac$est[,2,i]
    
    acf(my_residuals, 
        na.action = na.pass, main = paste("Monitoring Site ", i, ": PM10 Residual Autocorrelation", sep=""),
        lag.max = 120) 
  }
}
dev.off()

pdf("SO2_ResidualAutocorrelation.pdf", width = 12, height = 8)
for(i in 1:12){
  if(sum(is.na(air_tensor_std[,3,i])) == 0){
    my_residuals <-air_tensor_std[,3,i] - res_ac$est[,3,i]
    
    acf(my_residuals, 
        na.action = na.pass, main = paste("Monitoring Site ", i, ": SO2 Residual Autocorrelation", sep=""),
        lag.max = 120) 
  }
}
dev.off()

pdf("NO2_ResidualAutocorrelation.pdf", width = 12, height = 8)
for(i in 1:12){
  if(sum(is.na(air_tensor_std[,4,i])) == 0){
    my_residuals <-air_tensor_std[,4,i] - res_ac$est[,4,i]
    
    acf(my_residuals, 
        na.action = na.pass, main = paste("Monitoring Site ", i, ": NO2 Residual Autocorrelation", sep=""),
        lag.max = 120) 
  }
}
dev.off()

pdf("CO_ResidualAutocorrelation.pdf", width = 12, height = 8)
for(i in 1:12){
  if(sum(is.na(air_tensor_std[,5,i])) == 0){
    my_residuals <-air_tensor_std[,5,i] - res_ac$est[,5,i]
    
    acf(my_residuals, 
        na.action = na.pass, main = paste("Monitoring Site ", i, ": CO Residual Autocorrelation", sep=""),
        lag.max = 120) 
  }
}
dev.off()

pdf("O3_ResidualAutocorrelation.pdf", width = 12, height = 8)
for(i in 1:12){
  if(sum(is.na(air_tensor_std[,6,i])) == 0){
    my_residuals <-air_tensor_std[,6,i] - res_ac$est[,6,i]
    
    acf(my_residuals, 
        na.action = na.pass, main = paste("Monitoring Site ", i, ": O3 Residual Autocorrelation", sep=""),
        lag.max = 120) 
  }
}
dev.off()


## Run h-block cross validation with 6 folds (2-months validation data) and 7-day buffer before and after the validation data
## A rough search first
kcv_res <- kcv_hblock(tnsr=air_tensor_std, rank_grid=as.matrix(expand.grid(r1<-seq(3,6,by=1), r2<-seq(2,5, by=1))), 
               lambda_seq=c(0, 10^(1:6)), k=6, h=7, L0=NULL, D=D2, tol=0.1, max_iter=500, init=0)

kcv_res 
# 5    5 1e+05

# $MSE_mat
# [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
# [1,] 0.9997723 0.9588923 0.9847406 0.7184769 0.6407791 0.6187372 0.6127797
# [2,] 0.9997723 0.9585836 0.9833316 0.7192645 0.6410632 0.6188413 0.6127149
# [3,] 0.9997723 0.9585005 0.9836756 0.7196969 0.6404141 0.6178707 0.6120936
# [4,] 0.9997723 0.9588054 0.9840230 0.7204108 0.6405570 0.6179690 0.6121005
# [5,] 0.9997723 0.9583943 0.9550731 0.6662116 0.5970302 0.5836202 0.5803709
# [6,] 0.9997723 0.9593092 0.9716391 0.7054795 0.6106659 0.5816090 0.5762706
# [7,] 0.9997723 0.9588119 0.9709855 0.7036064 0.6082141 0.5785124 0.5744857
# [8,] 0.9997723 0.9588197 0.9721144 0.7046370 0.6086724 0.5789678 0.5746104
# [9,] 0.9997723 0.9614685 0.9905824 0.6810069 0.6241076 0.6119762 0.6088714
# [10,] 0.9997723 0.9586947 0.9811626 0.7175831 0.6162910 0.5843944 0.5884107
# [11,] 0.9997723 0.9582046 0.9818953 0.6936068 0.5877356 0.5591895 0.5759157
# [12,] 0.9997723 0.9582946 0.9818680 0.7060142 0.5878200 0.5597738 0.5760603
# [13,] 0.9997723 0.9615468 0.9913712 0.6877499 0.6307154 0.6185137 0.6145122
# [14,] 0.9997723 0.9587332 0.9838344 0.7359750 0.6432336 0.6062593 0.6030530
# [15,] 0.9997723 0.9572774 0.9807189 0.6931596 0.5720556 0.5489214 0.5781851
# [16,] 0.9997723 0.9571719 0.9811444 0.6955861 0.5770976 0.5496806 0.5782090
# 
# $SE_mat
# [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
# [1,] 0.1471379 0.1199308 0.1664337 0.1328582 0.1359453 0.1296471 0.1191261
# [2,] 0.1471379 0.1197760 0.1673987 0.1329126 0.1362189 0.1294794 0.1190410
# [3,] 0.1471379 0.1198110 0.1678503 0.1333496 0.1361874 0.1294096 0.1189952
# [4,] 0.1471379 0.1198818 0.1681698 0.1338589 0.1361421 0.1293786 0.1189910
# [5,] 0.1471379 0.1195050 0.1514941 0.1098717 0.1285847 0.1270739 0.1193131
# [6,] 0.1471379 0.1204276 0.1647501 0.1332993 0.1343839 0.1294026 0.1198216
# [7,] 0.1471379 0.1202419 0.1647671 0.1330014 0.1343216 0.1290112 0.1195455
# [8,] 0.1471379 0.1202656 0.1655185 0.1334834 0.1345483 0.1291848 0.1195874
# [9,] 0.1471379 0.1220458 0.1703881 0.1117483 0.1307000 0.1293736 0.1211691
# [10,] 0.1471379 0.1203226 0.1757937 0.1437883 0.1320897 0.1294153 0.1154160
# [11,] 0.1471379 0.1200193 0.1785618 0.1386551 0.1378977 0.1308346 0.1160718
# [12,] 0.1471379 0.1202338 0.1789374 0.1480165 0.1380300 0.1311971 0.1161612
# [13,] 0.1471379 0.1220575 0.1707257 0.1105061 0.1289009 0.1273606 0.1195055
# [14,] 0.1471379 0.1205578 0.1753605 0.1433107 0.1335128 0.1331603 0.1207566
# [15,] 0.1471379 0.1190214 0.1780226 0.1481879 0.1323908 0.1288621 0.1192167
# [16,] 0.1471379 0.1191379 0.1781612 0.1483990 0.1340736 0.1286786 0.1192181


## When r1=5. r2=5, lambda=1e+05 look into the variability explained by each component 
res_opt <-  mglram(tnsr = air_tensor_std, ranks = c(5, 5), init=0, D = D2,
                   lambda = 100000, max_iter = 500, tol = 1e-5, L0 = NULL)
res_opt$conv

tilde <- MakeIdent(res_opt$L, res_opt$G, res_opt$R)
L_tilde <- tilde$L_tilde
G_tilde <- tilde$G_tilde
R_tilde <- tilde$R_tilde

cum_var_L <- rep(NA, 5)
for (i in 1:5){
  comp <- array(NA, dim(air_tensor_std))
  for (j in 1:12){
    comp[ , , j] <- matrix(L_tilde[, 1:i], ncol=i) %*% matrix(G_tilde[1:i, , j], nrow=i) %*% t(R_tilde)
  }
  cum_var_L[i] <- sum(comp^2)/sum(air_tensor_std^2)
}

cum_var_L # 0.3530622 0.4417132 0.4918228 0.5096191 0.5184378

cum_var_L0 <- c(0, cum_var_L[1:4])
sep_var_L <- cum_var_L - cum_var_L0
sep_var_L # 0.353062161 0.088651020 0.050109643 0.017796238 0.008818764

cum_var_R <- rep(NA, 5)
for (i in 1:5){
  comp <- array(NA, dim(air_tensor_std))
  for (j in 1:12){
    comp[ , , j] <- L_tilde %*% matrix(G_tilde[,1:i, j], ncol=i) %*% matrix(t(R_tilde[, 1:i]), nrow=i)
  }
  cum_var_R[i] <- sum(comp^2)/sum(air_tensor_std^2)
}

cum_var_R # 0.4049534 0.4558021 0.4926196 0.5089911 0.5184378

cum_var_R0 <- c(0, cum_var_R[1:4])
sep_var_R <- cum_var_R - cum_var_R0

sep_var_R # 0.404953398 0.050848710 0.036817456 0.016371584 0.009446679

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

#ggsave("./air_quality_var_L.pdf", dpi=600, width=6, height=4)

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

#ggsave("./air_quality_var_R.pdf", dpi=600, width=6, height=4)


kcv_fine <- kcv_hblock(tnsr=air_tensor_std, rank_grid=as.matrix(expand.grid(r1<-c(3), r2<-c(3))),
                  lambda_seq=c(80000, 90000, 100000,200000, 300000, 400000, 500000), k=6, h=7, L0=NULL, D=D2, tol=0.1, max_iter=500, init=0)

# $MSE_mat
# [,1]      [,2]      [,3]     [,4]      [,5]      [,6]      [,7]
# [1,] 0.5859994 0.5847279 0.5836267 0.577378 0.5755042 0.5750465 0.5753571
# 
# $SE_mat
# [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
# [1,] 0.1279563 0.1275007 0.1270707 0.1243994 0.1230222 0.1220771 0.1213553
# 
# $opt_para
# [,1] [,2]  [,3]
# [1,]    3    3 4e+05

res <- mglram(tnsr = air_tensor_std, ranks = c(3, 3), init=0, D = D2,
              lambda = 400000, max_iter = 500, tol = 1e-5, L0 = NULL)

tilde <- MakeIdent(L=res$L, G=res$G, R=res$R)
L_tilde <- tilde$L_tilde
R_tilde <- tilde$R_tilde
G_tilde <- tilde$G_tilde

sum(res$est^2)/sum(air_tensor_std^2) #  0.4610038

month_mids  <- c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349)
#month_labs  <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
month_labs <- c("1","2","3","4","5","6","7","8","9","10","11","12")

dataL = data.frame(loading = c(L_tilde[, 1], L_tilde[, 2], L_tilde[, 3]),
                   component = c(rep("1st", 366), rep("2nd", 366), rep("3rd", 366)),
                   t = rep(1:366, 3))

pL = dataL %>%
  ggplot(aes(x = t, y = loading)) +
  geom_vline(xintercept = month_mids, color = "grey80", linewidth = 0.4) +  # month midlines
  geom_line(alpha = 0.5) +
  geom_point(size = 0.1) +
  geom_hline(yintercept = 0, color = "red") +
  scale_x_continuous(breaks = month_mids, labels = month_labs) +
  xlab("Month") +
  facet_grid(~ component) +
  ylab("Coefficient") +
  scale_y_continuous(limits = c(-0.1, 0.1)) +
  theme_bw() +
  theme(
    axis.title = element_text(hjust=0.5, size=15),
    axis.text.y = element_text(hjust=0.5, size=12),
    axis.text.x = element_text(hjust=0.5, size=12),
    strip.text = element_text(size = 15),
    legend.text = element_text(size = 15),       
    legend.title = element_text(size = 15),
    strip.placement = "outside",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(color = NA)
  )
pL

#ggsave("./air_quality_L.pdf", width=8, height=3)

## Plot effect of time components on different variables
n = 12

Ri = array(NA, dim = c(3, 6, n))
for (i in 1:n){
  Ri[, , i] = G_tilde[, , i] %*% t(R_tilde)
}

RiL1 = Ri[1, , ] 
apply(RiL1, 1, sd)  # 0.7813149 0.9683754 2.6105632 1.2363402 1.1976043 1.1010215

RiL2 = Ri[2, , ]
apply(RiL2, 1, sd)  # 3.405264  4.512273  3.816563 11.528460  4.106275  2.741775

RiL3 = Ri[3, , ]
apply(RiL3, 1, sd)  # 0.4511773 0.5815731 1.2994171 1.8084842 0.5735618 0.4721721

## PM2.5: Mean=79.71, SD=35.73
PM25mean = rep(79.71, 366) + (L_tilde[, 1] * mean(RiL1[1,]) + L_tilde[, 2] * mean(RiL2[1,]) + L_tilde[, 3] * mean(RiL3[1,]))*35.73

dataLeffects = data.frame(PM25mean = rep(PM25mean, 6), 
                          value = c(PM25mean + L_tilde[, 1] *0.78*35.73, PM25mean - L_tilde[, 1] * 0.78*35.73, 
                                    PM25mean + L_tilde[, 2] *3.41*35.73, PM25mean - L_tilde[, 2] * 3.41*35.73,
                                    PM25mean + L_tilde[, 3] *0.45*35.73, PM25mean - L_tilde[, 3] * 0.45*35.73),
                          component = c(rep("1st", 366 * 2), rep("2nd", 366 * 2), rep("3rd",366*2)),
                          sign = rep(rep(c("positive", "negative"), each = 366), 3),
                          t = rep(1:366, 6))

pdataL = dataLeffects %>%
  ggplot(aes(x = t, y = value, group = sign, col = sign)) +
  geom_vline(xintercept = month_mids, color = "grey80", linewidth = 0.2) +
  geom_line(alpha = 2) +
  xlab("Month") + 
  facet_grid(~ component) + 
  ylab("PM2.5") + geom_line(aes(x = t, y = PM25mean), col = "black") +
  scale_x_continuous(breaks = month_mids, labels = month_labs) +
  theme_bw() +
  theme(
    axis.title = element_text(hjust=0.5, size=15),
    axis.text.y = element_text(hjust=0.5, size=12),
    axis.text.x = element_text(hjust=0.5, size=12),
    strip.text = element_text(size = 15),
    legend.text = element_text(size = 15),       
    legend.title = element_text(size = 15),
    strip.placement = "outside",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(color = NA),
    legend.position = "bottom"
  )

pdataL

#ggsave("./PM25.pdf", width=8, height=3)

## PM10: Mean=104.52, SD=40.59 
PM10mean = rep(104.52, 366) + (L_tilde[, 1] * mean(RiL1[2,]) + L_tilde[, 2] * mean(RiL2[2,]) + L_tilde[, 3] * mean(RiL3[2,]))*40.59

dataLeffects = data.frame(PM10mean = rep(PM10mean, 6), 
                          value = c(PM10mean + L_tilde[, 1] * 0.97*40.59, PM10mean - L_tilde[, 1] * 0.97*40.59, 
                                    PM10mean + L_tilde[, 2] *4.51*40.59, PM10mean - L_tilde[, 2] * 4.51*40.59,
                                    PM10mean + L_tilde[, 3] *0.58*40.59, PM10mean - L_tilde[, 3] * 0.58*40.59),
                          component = c(rep("1st", 366 * 2), rep("2nd", 366 * 2), rep("3rd",366*2)),
                          sign = rep(rep(c("positive", "negative"), each = 366), 3),
                          t = rep(1:366, 6))

pdataL = dataLeffects %>%
  ggplot(aes(x = t, y = value, group = sign, col = sign)) +
  geom_vline(xintercept = month_mids, color = "grey80", linewidth = 0.2) +
  geom_line(alpha = 2) +
  xlab("Month") + 
  facet_grid(~ component) + 
  ylab("PM10") + geom_line(aes(x = t, y = PM10mean), col = "black") +
  scale_x_continuous(breaks = month_mids, labels = month_labs) +
  theme_bw() +
  theme(
    axis.title = element_text(hjust=0.5, size=15),
    axis.text.y = element_text(hjust=0.5, size=12),
    axis.text.x = element_text(hjust=0.5, size=12),
    strip.text = element_text(size = 15),
    legend.text = element_text(size = 15),       
    legend.title = element_text(size = 15),
    strip.placement = "outside",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(color = NA),
    legend.position = "bottom"
  )

pdataL

#ggsave("./PM10.pdf", width=8, height=3)

## SO2: Mean=15.80, SD=11.84
SO2mean = rep(15.8, 366) + (L_tilde[, 1] * mean(RiL1[3,]) + L_tilde[, 2] * mean(RiL2[3,]) + L_tilde[, 3] * mean(RiL3[3,]))*11.84

dataLeffects = data.frame(SO2mean = rep(SO2mean, 6), 
                          value = c(SO2mean + L_tilde[, 1] * 2.61*11.84, SO2mean - L_tilde[, 1] * 2.61*11.84, 
                                    SO2mean + L_tilde[, 2] *3.82*11.84, SO2mean - L_tilde[, 2] * 3.82*11.84,
                                    SO2mean + L_tilde[, 3] *1.30*11.84, SO2mean - L_tilde[, 3] * 1.30*11.84),
                          component = c(rep("1st", 366 * 2), rep("2nd", 366 * 2), rep("3rd",366*2)),
                          sign = rep(rep(c("positive", "negative"), each = 366), 3),
                          t = rep(1:366, 6))

pdataL = dataLeffects %>%
  ggplot(aes(x = t, y = value, group = sign, col = sign)) +
  geom_vline(xintercept = month_mids, color = "grey80", linewidth = 0.2) +
  geom_line(alpha = 2) +
  xlab("Month") + 
  facet_grid(~ component) + 
  ylab("SO2") + geom_line(aes(x = t, y = SO2mean), col = "black") +
  scale_x_continuous(breaks = month_mids, labels = month_labs) +
  theme_bw() +
  theme(
    axis.title = element_text(hjust=0.5, size=15),
    axis.text.y = element_text(hjust=0.5, size=12),
    axis.text.x = element_text(hjust=0.5, size=12),
    strip.text = element_text(size = 15),
    legend.text = element_text(size = 15),       
    legend.title = element_text(size = 15),
    strip.placement = "outside",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(color = NA),
    legend.position = "bottom"
  )

pdataL

#ggsave("./SO2.pdf", width=8, height=3)

## NO2: Mean=50.55, SD=18.91
NO2mean = rep(50.55, 366) + (L_tilde[, 1] * mean(RiL1[4,]) + L_tilde[, 2] * mean(RiL2[4,]) + L_tilde[, 3] * mean(RiL3[4,]))*18.91

dataLeffects = data.frame(NO2mean = rep(NO2mean, 6), 
                          value = c(NO2mean + L_tilde[, 1] * 1.24*18.91, NO2mean - L_tilde[, 1] * 1.24*18.91, 
                                    NO2mean + L_tilde[, 2] *11.53*18.91, NO2mean - L_tilde[, 2] * 11.53*18.91,
                                    NO2mean + L_tilde[, 3] *1.81*18.91, NO2mean - L_tilde[, 3] * 1.81 *18.91),
                          component = c(rep("1st", 366 * 2), rep("2nd", 366 * 2), rep("3rd",366*2)),
                          sign = rep(rep(c("positive", "negative"), each = 366), 3),
                          t = rep(1:366, 6))

pdataL = dataLeffects %>%
  ggplot(aes(x = t, y = value, group = sign, col = sign)) +
  geom_vline(xintercept = month_mids, color = "grey80", linewidth = 0.2) +
  geom_line(alpha = 2) +
  xlab("Month") + 
  facet_grid(~ component) + 
  ylab("NO2") + geom_line(aes(x = t, y = NO2mean), col = "black") +
  scale_x_continuous(breaks = month_mids, labels = month_labs) +
  theme_bw() +
  theme(
    axis.title = element_text(hjust=0.5, size=15),
    axis.text.y = element_text(hjust=0.5, size=12),
    axis.text.x = element_text(hjust=0.5, size=12),
    strip.text = element_text(size = 15),
    legend.text = element_text(size = 15),       
    legend.title = element_text(size = 15),
    strip.placement = "outside",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(color = NA),
    legend.position = "bottom"
  )

pdataL

#ggsave("./NO2.pdf", width=8, height=3)

## CO: Mean=1232.44, SD=633.19
COmean = rep(1232.44, 366) + (L_tilde[, 1] * mean(RiL1[5,]) + L_tilde[, 2] * mean(RiL2[5,]) + L_tilde[, 3] * mean(RiL3[5,]))*633.19

dataLeffects = data.frame(COmean = rep(COmean, 6), 
                          value = c(COmean + L_tilde[, 1] * 1.20*633.19, COmean - L_tilde[, 1] * 1.20*633.19, 
                                    COmean + L_tilde[, 2] *4.11*633.19, COmean - L_tilde[, 2] * 4.11*633.19,
                                    COmean + L_tilde[, 3] *0.57*633.19, COmean - L_tilde[, 3] * 0.57*633.19),
                          component = c(rep("1st", 366 * 2), rep("2nd", 366 * 2), rep("3rd",366*2)),
                          sign = rep(rep(c("positive", "negative"), each = 366), 3),
                          t = rep(1:366, 6))

pdataL = dataLeffects %>%
  ggplot(aes(x = t, y = value, group = sign, col = sign)) +
  geom_vline(xintercept = month_mids, color = "grey80", linewidth = 0.2) +
  geom_line(alpha = 2) +
  xlab("Month") + 
  facet_grid(~ component) + 
  ylab("CO") + geom_line(aes(x = t, y = COmean), col = "black") +
  scale_x_continuous(breaks = month_mids, labels = month_labs) +
  theme_bw() +
  theme(
    axis.title = element_text(hjust=0.5, size=15),
    axis.text.y = element_text(hjust=0.5, size=12),
    axis.text.x = element_text(hjust=0.5, size=12),
    strip.text = element_text(size = 15),
    legend.text = element_text(size = 15),       
    legend.title = element_text(size = 15),
    strip.placement = "outside",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(color = NA),
    legend.position = "bottom"
  )

pdataL

#ggsave("./CO.pdf", width=8, height=3)

## O3
O3mean = rep(57.18, 366) + (L_tilde[, 1] * mean(RiL1[6,]) + L_tilde[, 2] * mean(RiL2[6,]) + L_tilde[, 3] * mean(RiL3[6,]))*31.83 

dataLeffects = data.frame(O3mean = rep(O3mean, 6), 
                          value = c(O3mean + L_tilde[, 1] * 1.10*31.83, O3mean - L_tilde[, 1] * 1.10*31.83, 
                                    O3mean + L_tilde[, 2] *2.74*31.83, O3mean - L_tilde[, 2] * 2.74*31.83,
                                    O3mean + L_tilde[, 3] *0.47*31.83, O3mean - L_tilde[, 3] * 0.47*31.83),
                          component = c(rep("1st", 366 * 2), rep("2nd", 366 * 2), rep("3rd",366*2)),
                          sign = rep(rep(c("positive", "negative"), each = 366), 3),
                          t = rep(1:366, 6))

pdataL = dataLeffects %>%
  ggplot(aes(x = t, y = value, group = sign, col = sign)) +
  geom_vline(xintercept = month_mids, color = "grey80", linewidth = 0.2) +
  geom_line(alpha = 2) +
  xlab("Month") + 
  facet_grid(~ component) + 
  ylab("O3") + geom_line(aes(x = t, y = O3mean), col = "black") +
  scale_x_continuous(breaks = month_mids, labels = month_labs) +
  theme_bw() +
  theme(
    axis.title = element_text(hjust=0.5, size=15),
    axis.text.y = element_text(hjust=0.5, size=12),
    axis.text.x = element_text(hjust=0.5, size=12),
    strip.text = element_text(size = 15),
    legend.text = element_text(size = 15),       
    legend.title = element_text(size = 15),
    strip.placement = "outside",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(color = NA),
    legend.position = "bottom"
  )

pdataL

#ggsave("./O3.pdf", width=8, height=3)
