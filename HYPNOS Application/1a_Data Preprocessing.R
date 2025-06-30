load("Application/filtered_ABPM_visit1_p.Rda")

library(ggplot2)
library(bp)
library(dplyr)

plots = bp_ts_plots(filtered_abpm, subj = 70008)
plots # Blood pressure by Date/Time and Blood pressure by hours (may have repeated measurements at one time)

# Check outliers in DBP, SBP, and MAP 
# When SBP/DBP are abnormal, the MAP will also be abnormal since it is like an average of SBP and DBP. Therefore, only need to check SBP and DBP
# According to process_data in bp package, we want to set to NA the ones with SBP above 240 or below 50 (SUL and SLL)
filtered_abpm[filtered_abpm$SBP > 240, ] #none
filtered_abpm[filtered_abpm$SBP < 50, ] #none 
# According to process_data in bp package, we want to set to NA the ones with DBP above 140 or below 40 (DUL and DLL)
filtered_abpm[filtered_abpm$DBP > 140, ] # 2 measurements, one for 70179 and one for 70638
filtered_abpm[filtered_abpm$DBP < 40, ] # many measurements

# Check 70179
plots = bp_ts_plots(filtered_abpm, subj = 70179)
plots[[1]] # There is an outlier for both SBP and DBP
# Check 70638 
plots = bp_ts_plots(filtered_abpm, subj = 70638)
plots[[1]] # There is an outlier for DBP only

# Check IDs with DBP < 40
plots = bp_ts_plots(filtered_abpm, subj = 70080)
plots[[1]] 
plots = bp_ts_plots(filtered_abpm, subj = 70107)
plots[[1]]
plots = bp_ts_plots(filtered_abpm, subj = 70109)
plots[[1]] 
plots = bp_ts_plots(filtered_abpm, subj = 70153)
plots[[1]]
plots = bp_ts_plots(filtered_abpm, subj = 70362)
plots[[1]]
plots = bp_ts_plots(filtered_abpm, subj = 70424)
plots[[1]] 
plots = bp_ts_plots(filtered_abpm, subj = 70488)
plots[[1]]
plots = bp_ts_plots(filtered_abpm, subj = 70545)
plots[[1]] 
plots = bp_ts_plots(filtered_abpm, subj = 70574)
plots[[1]] 
plots = bp_ts_plots(filtered_abpm, subj = 70694)
plots[[1]] # SBP and DBP both abnormal (too low)
plots = bp_ts_plots(filtered_abpm, subj = 90034)
plots[[1]] # SBP and DBP both abnormal (too low)
plots = bp_ts_plots(filtered_abpm, subj = 90070)
plots[[1]] 
plots = bp_ts_plots(filtered_abpm, subj = 90082)
plots[[1]]
plots = bp_ts_plots(filtered_abpm, subj = 90099)
plots[[1]]

# Check outliers in HR
# According to process_data in bp package, the HR should be below 220 and above 27
filtered_abpm[filtered_abpm$HR > 220, ] # none
filtered_abpm[filtered_abpm$HR < 27, ] #none 

# Remove outliers in SBP, DBP, MAP, HR (no outlier)
filtered_abpm_new =  filtered_abpm %>%
  mutate(SBP = ifelse(((SBP > 200)&(ID == 70179))| ((DBP < 40)&(ID == 70694)), NA, SBP),
         DBP = ifelse(DBP > 140 | DBP < 40, NA, DBP), 
         MAP = ifelse(is.na(SBP) | is.na(DBP), NA, MAP))

# save(filtered_abpm_new, file = "filtered_ABPM_new.Rda")

hourly_abpm  = filtered_abpm_new %>%
  group_by(ID, HOUR) %>%
  summarize(SBPh = mean(SBP, na.rm = TRUE), DBPh = mean(DBP, na.rm = TRUE), MAPh = mean(MAP, na.rm = TRUE), HRh = mean(HR, na.rm=TRUE))

# save(hourly_abpm, file = "hourly_ABPM_v1.Rda")


