load("Application/hourly_ABPM_v1.Rda")

library(rTensor)
library(dplyr)
library(tidyr)
library(refund)

subjects = unique(hourly_abpm$ID)
n = length(subjects)

## Assign SBP, DBP, MAP, HR into matrices 
library(tidyr)
SBP_wide = pivot_wider(hourly_abpm, id_cols = "ID", names_from = "HOUR", values_from = "SBPh")
DBP_wide = pivot_wider(hourly_abpm, id_cols = "ID", names_from = "HOUR", values_from = "DBPh")
MAP_wide = pivot_wider(hourly_abpm, id_cols = "ID", names_from = "HOUR", values_from = "MAPh")
HR_wide = pivot_wider(hourly_abpm, id_cols = "ID", names_from = "HOUR", values_from = "HRh")

## Remove ID column
SBPm = SBP_wide[, -1]
DBPm = DBP_wide[, -1]
MAPm = MAP_wide[, -1]
HRm = HR_wide[, -1]

## Rearrange columns corresponding to hours
## First order from 0 to 23
new_order = order(as.numeric(names(SBP_wide)[-1]))
SBPm = SBPm[, new_order]
DBPm = DBPm[, new_order]
MAPm = MAPm[, new_order]
HRm = HRm[, new_order]

## Next put the 0 to 11 in the end
SBPm = SBPm[, c(13:24, 1:12)]
DBPm = DBPm[, c(13:24, 1:12)]
MAPm = MAPm[, c(13:24, 1:12)]
HRm = HRm[, c(13:24, 1:12)]

## Check order
names(SBPm)
names(DBPm)
names(MAPm)
names(HRm)

SBPnew = as.matrix(SBPm)
DBPnew = as.matrix(DBPm)
MAPnew = as.matrix(MAPm)
HRnew = as.matrix(HRm)

mean(SBPnew,na.rm=T) # 128.7944
mean(DBPnew,na.rm=T) # 70.58507
mean(MAPnew,na.rm=T) # 90.33847
mean(HRnew,na.rm=T) # 76.23369

sd(DBPnew, na.rm = T) # 12.32371 
sd(SBPnew, na.rm = T) # 18.30611
sd(MAPnew, na.rm = T) # 13.49582
sd(HRnew, na.rm = T) # 12.85726

## Put hourly data into a tensor format 
## n subjects by 24 hours of the day by 4 measurements (SBP/DBP/MAP/HR)
BP_visit1 <- array(NA, dim = c(n, 24, 4))

## Normalize the data
BP_visit1[, , 1] = (DBPnew - mean(DBPnew,na.rm=T))/sd(DBPnew, na.rm = T) 
BP_visit1[, , 2] = (SBPnew - mean(SBPnew,na.rm=T))/sd(SBPnew, na.rm = T) 
BP_visit1[, , 3] = (MAPnew - mean(MAPnew,na.rm=T))/sd(MAPnew, na.rm = T) 
BP_visit1[, , 4] = (HRnew - mean(HRnew,na.rm=T))/sd(HRnew, na.rm = T)
BP_visit1 = as.tensor(BP_visit1)

newBP = array(NA, dim = c(24, 4, 207))
newBP[, 1, ] = t(BP_visit1@data[ , , 1])
newBP[, 2, ] = t(BP_visit1@data[ , , 2])
newBP[, 3, ] = t(BP_visit1@data[ , , 3])
newBP[, 4, ] = t(BP_visit1@data[ , , 4])
newBP = as.tensor(newBP)

sum(is.na(newBP@data))

missing_normalized_4 <- newBP

# save(missing_normalized_4,file="missing_normalized_4.Rda")

newBP_3 = array(NA, dim = c(24, 3, 207))
newBP_3[, 1, ] = t(BP_visit1@data[ , , 1])
newBP_3[, 2, ] = t(BP_visit1@data[ , , 2])
newBP_3[, 3, ] = t(BP_visit1@data[ , , 4])
newBP_3 = as.tensor(newBP_3)

## Number and proportion of missing elements
sum(is.na(newBP_3@data)) 
sum(is.na(newBP_3@data))/(24*3*207) # 15.9%

## Missing rate for each patient
missing_rate_vec <- rep(0, 207)
for (i in 1:207){
  missing_rate_vec[i] <- sum(is.na(newBP_3@data[ , , i]))/(24*3)
}
sum(missing_rate_vec==0)/207 # 18.4%
min(missing_rate_vec[missing_rate_vec!=0]) # 4.2%
max(missing_rate_vec[missing_rate_vec!=0]) # 79.2%

missing_normalized_3 <- newBP_3

# save(missing_normalized_3,file="missing_normalized_3.Rda")

## Smooth using fPCA to deal with NAs
## SBP
outfpcaSBP <- fpca.sc(Y = as.matrix(SBPm), pve = 0.9) 
outfpcaSBP$evalues

## DBP
outfpcaDBP <- fpca.sc(Y = as.matrix(DBPm), pve = 0.9) 
outfpcaDBP$evalues

## MAP 
outfpcaMAP <- fpca.sc(Y = as.matrix(MAPm), pve = 0.9) 
outfpcaMAP$evalues

## HR
outfpcaHR <- fpca.sc(Y = as.matrix(HRm), pve = 0.9) 
outfpcaHR$evalues

SBPnew = as.matrix(SBPm)
SBPnew[is.na(SBPnew)] = outfpcaSBP$Yhat[is.na(SBPnew)] 
DBPnew = as.matrix(DBPm)
DBPnew[is.na(DBPnew)] = outfpcaDBP$Yhat[is.na(DBPnew)] 
MAPnew = as.matrix(MAPm)
MAPnew[is.na(MAPnew)] = outfpcaMAP$Yhat[is.na(MAPnew)] 
HRnew = as.matrix(HRm)
HRnew[is.na(HRnew)] = outfpcaHR$Yhat[is.na(HRnew)] 

## Put hourly data into a tensor format 
## n subjects by 24 hours of the day by 4 measurements (SBP/DBP/MAP/HR)
BP_visit1_4 <- array(NA, dim = c(n, 24, 4))
## n subjects by 24 hours of the day by 3 measurements (SBP/DBP/HR)
BP_visit1_3 <- array(NA, dim = c(n, 24, 3))

## Overall mean across all times
mean_original_DBP <- mean(as.matrix(DBPm),na.rm = T)
mean_original_SBP <- mean(as.matrix(SBPm),na.rm = T)
mean_original_MAP <- mean(as.matrix(MAPm),na.rm = T)
mean_original_HR <- mean(as.matrix(HRm),na.rm = T)

sd_original_DBP <- sd(as.matrix(DBPm),na.rm = T)
sd_original_SBP <-sd(as.matrix(SBPm),na.rm = T)
sd_original_MAP <- sd(as.matrix(MAPm),na.rm = T)
sd_original_HR <-sd(as.matrix(HRm),na.rm = T)

BP_visit1_4[, , 1] = (DBPnew - mean_original_DBP)/sd_original_DBP 
BP_visit1_4[, , 2] = (SBPnew - mean_original_SBP)/sd_original_SBP 
BP_visit1_4[, , 3] = (MAPnew - mean_original_MAP)/sd_original_MAP 
BP_visit1_4[, , 4] = (HRnew - mean_original_HR)/sd_original_HR 
BP_visit1_4 = as.tensor(BP_visit1_4)

newBP_filled = array(NA, dim = c(24, 4, 207))
newBP_filled[, 1, ] = t(BP_visit1_4@data[ , , 1])
newBP_filled[, 2, ] = t(BP_visit1_4@data[ , , 2])
newBP_filled[, 3, ] = t(BP_visit1_4@data[ , , 3])
newBP_filled[, 4, ] = t(BP_visit1_4@data[ , , 4])
newBP_filled = as.tensor(newBP_filled)

fpca_filled <- newBP_filled
# save(fpca_filled,file="fpca_filled.Rda")


