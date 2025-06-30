load("./Application/missing_normalized_3.Rda")
load("./Application/fpca_filled.Rda") 
load("./Application/AlgorithmResIdent.Rda") # L_tilde, R_tilde, G_tilde

library(rTensor)
library(dplyr)
library(ggplot2)
library(car)
library(tidyverse)
library(broom)
library(readxl)

# Plot coefficient for L components
dataL = data.frame(loading = c(L_tilde[, 1], L_tilde[, 2], L_tilde[, 3]),
                   component = c(rep("1st", 24), rep("2nd", 24), rep("3rd", 24)),
                   hour = rep(12:35, 3))

pL = dataL %>%
  ggplot(aes(x = hour, y = loading)) + geom_line(alpha = 0.5) + geom_point() + geom_hline(yintercept = 0, color = "red") +
  scale_x_continuous(breaks = c(12, 15, 18, 21, 24, 27, 30, 33, 36), labels = c(12, 15, 18, 21, 0, 3, 6, 9, 12)) +
  xlab("Hour of the day") + 
  facet_grid(~ component) + 
  ylab("Coefficient") +
  theme(text = element_text(size = 15)) +
  scale_y_continuous(limits = c(-0.4, 0.4))

pL

# ggsave("./Application/L_comp.pdf", dpi=600, width=8, height=2.5)

# Plot effect of time components 
n = 207

Ri = array(NA, dim = c(3, 3, n))
for (i in 1:n){
  Ri[, , i] = G_tilde[, , i] %*% t(R_tilde)
}

RiL1 = Ri[1, , ] 
apply(RiL1, 1, sd) # 2.789581 2.997533 3.661661 (sd for DBP, SBP, and HR for 1st time component)

RiL2 = Ri[2, , ]
apply(RiL2, 1, sd) # 1.379921 1.401255 1.444989 (sd for DBP, SBP, and HR for 2nd time component)

RiL3 = Ri[3, , ]
apply(RiL3, 1, sd) # 1.067016 1.115904 1.072485 (sd for DBP, SBP, and HR for 3rd time component)


DBPmean = rep(70, 24) + (L_tilde[, 1] * mean(RiL1[1,]) + L_tilde[, 2] * mean(RiL2[1,]) + L_tilde[, 3] * mean(RiL3[1,])) * 12 # 70 is the mean of DBP, 12 is the sd of DBP, doing this because we normalized the data initially
SBPmean = rep(129, 24) + (L_tilde[, 1] * mean(RiL1[2,]) + L_tilde[, 2] * mean(RiL2[2,]) + L_tilde[, 3] * mean(RiL3[2,])) * 18 # 129 is the mean of SBP, 18 is the sd of SBP
HRmean = rep(76, 24) + (L_tilde[, 1] * mean(RiL1[3,]) + L_tilde[, 2] * mean(RiL2[3,]) + L_tilde[, 3] * mean(RiL3[3,])) * 13 # 76 is the mean of HR, 13 is the sd of HR 

# Plot effect of time components on DBP 
dataLeffects = data.frame(DBPmean = rep(DBPmean, 6), 
                          value = c(DBPmean + L_tilde[, 1] * 2.8*12, DBPmean - L_tilde[, 1] * 2.8*12, 
                                    DBPmean + L_tilde[, 2] *1.4*12, DBPmean - L_tilde[, 2] * 1.4*12, 
                                    DBPmean + L_tilde[, 3] * 1.1*12, DBPmean - L_tilde[, 3] * 1.1*12),
                          component = c(rep("1st", 24 * 2), rep("2nd", 24 * 2), rep("3rd", 24 * 2)),
                          sign = rep(rep(c("positive", "negative"), each = 24), 3),
                          hour = rep(12:35, 6))

pdataL = dataLeffects %>%
  ggplot(aes(x = hour, y = value, group = sign, col = sign)) + geom_line(alpha = 0.5) +
  scale_x_continuous(breaks = c(12, 15, 18, 21, 24, 27, 30, 33, 36), labels = c(12, 15, 18, 21, 0, 3, 6, 9, 12)) +
  xlab("Hour of the day") + 
  facet_grid(~ component) + 
  ylab("DBP") + geom_line(aes(x = hour, y = DBPmean), col = "black") +
  theme(text = element_text(size = 15)) 

pdataL

# ggsave("./Application/DBP.pdf", dpi=600, width=8, height=2.5)


# Plot effect of time components on SBP 
dataLeffects = data.frame(SBPmean = rep(SBPmean, 6), 
                          value = c(SBPmean + L_tilde[, 1] * 3*18, SBPmean - L_tilde[, 1] * 3*18, 
                                    SBPmean + L_tilde[, 2] * 1.4*18, SBPmean - L_tilde[, 2] * 1.4*18, 
                                    SBPmean + L_tilde[, 3] * 1.1*18, SBPmean - L_tilde[, 3] * 1.1*18),
                          component = c(rep("1st", 24 * 2), rep("2nd", 24 * 2), rep("3rd", 24 * 2)),
                          sign = rep(rep(c("positive", "negative"), each = 24), 3),
                          hour = rep(12:35, 6))

pdataL = dataLeffects %>%
  ggplot(aes(x = hour, y = value, group = sign, col = sign)) + geom_line(alpha = 0.5) +
  scale_x_continuous(breaks = c(12, 15, 18, 21, 24, 27, 30, 33, 36), labels = c(12, 15, 18, 21, 0, 3, 6, 9, 12)) +
  xlab("Hour of the day") + facet_grid(~ component) + ylab("SBP") + geom_line(aes(x = hour, y = SBPmean), col = "black") +
  theme(text = element_text(size = 15)) 

pdataL

# ggsave("./Application/SBP.pdf", dpi=600,  width=8, height=2.5)

# Plot effect of time components on HR
dataLeffects = data.frame(HRmean = rep(HRmean, 6), 
                          value = c(HRmean + L_tilde[, 1] * 3.7*13, HRmean - L_tilde[, 1] * 3.7*13, 
                                    HRmean + L_tilde[, 2] * 1.4*13, HRmean - L_tilde[, 2] * 1.4*13, 
                                    HRmean + L_tilde[, 3] * 1.1*13, HRmean - L_tilde[, 3] * 1.1*13),
                          component = c(rep("1st", 24 * 2), rep("2nd", 24 * 2), rep("3rd", 24 * 2)),
                          sign = rep(rep(c("positive", "negative"), each = 24), 3),
                          hour = rep(12:35, 6))

pdataL = dataLeffects %>%
  ggplot(aes(x = hour, y = value, group = sign, col = sign)) + geom_line(alpha = 0.5) +
  scale_x_continuous(breaks = c(12, 15, 18, 21, 24, 27, 30, 33, 36), labels = c(12, 15, 18, 21, 0, 3, 6, 9, 12)) +
  xlab("Hour of the day") + facet_grid(~ component) + ylab("HR") + geom_line(aes(x = hour, y = HRmean), col = "black") +
  theme(text = element_text(size = 15)) 

pdataL

# ggsave("./Application/HR.pdf", dpi=600, width=8, height=2.5)

# What R looks like
R_tilde[, 1] # first component of R tilde
R_tilde[, 2] # second component of R tilde

## 37/100 DBP + 33/50 SBP + 29/100 HR (Seem to be joint effect of the three)
R_tilde[1, 1]/sum(R_tilde[, 1]) # DBP: 37/100 
R_tilde[2, 1]/sum(R_tilde[, 1]) # SBP: 33/100 
R_tilde[3, 1]/sum(R_tilde[, 1]) # HR: 30/100 = 3/10

## -13/100 DBP - 33/100 SBP + 27/50 HR (Contrast between blood pressure and heart rate?)
R_tilde[1, 2]/sum(abs(R_tilde[, 2])) # DBP: -13/100
R_tilde[2, 2]/sum(abs(R_tilde[, 2])) # SBP: -33/100 
R_tilde[3, 2]/sum(abs(R_tilde[, 2])) # HR: 54/100 = 27/50

# Plot Estimated Curves for DBP
est <- ttl(as.tensor(G_tilde), list_mat=list(L_tilde, R_tilde), ms=c(1,2))

for (i in 1:207){
  plot(missing_normalized_3@data[ , 1, i], main=paste(i,"th patient"), xlab="Hour", ylab="DBP", ylim=c(-1,1)) # original data
  lines(fpca_filled@data[ , 1, i], col="green") # data filled by fPCA
  lines(est@data[ , 1, i], col="blue") # estimated curve by our algorithm
}

# Plot Estimated Curves for SBP
for (i in 1:207){
  plot(missing_normalized_3@data[ , 2, i], main=paste(i,"th patient"), xlab="Hour", ylab="SBP", ylim=c(-2.5,2))
  lines(fpca_filled@data[ , 2, i], col="green")
  lines(est@data[ , 2, i], col="blue")
}

# Plot Estimated Curves for HR
for (i in 1:207){
  plot(missing_normalized_3@data[ , 3, i], main=paste(i,"th patient"), xlab="Hour", ylab="HR")
  lines(fpca_filled@data[ , 4, i], col="green")
  lines(est@data[ , 3, i], col="blue")
}
