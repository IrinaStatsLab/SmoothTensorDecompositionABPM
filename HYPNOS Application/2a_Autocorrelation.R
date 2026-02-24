load("./missing_normalized_3.Rda")

library(SmoothHOOI)

# Use parameters for moderate smoothing
D2 <- SecDiffMat(24)
res <- mglram(tnsr = missing_normalized_3@data, ranks = c(3, 2), init=0, D = D2,
              lambda = 50, max_iter = 500, tol = 1e-5, L0 = NULL)

# Check subjects with complete data
pdf("DBP_ResidualAutocorrelation.pdf", width = 12, height = 8)
for(i in 1:207){
  if(sum(is.na(missing_normalized_3@data[,1,i])) == 0){
    my_residuals <- missing_normalized_3@data[,1,i] - res$est[,1,i]
    
    acf(my_residuals, 
        na.action = na.pass, main = paste("Subject ", i, ": DBP Residual Autocorrelation", sep=""),
        lag.max = 24) 
  }
}
dev.off()

pdf("SBP_ResidualAutocorrelation.pdf", width = 12, height = 8)
for(i in 1:207){
  if(sum(is.na(missing_normalized_3@data[,2,i])) == 0){
    my_residuals <- missing_normalized_3@data[,2,i] - res$est[,2,i]
    
    acf(my_residuals, 
        na.action = na.pass, main = paste("Subject ", i, ": SBP Residual Autocorrelation", sep=""),
        lag.max = 24) 
  }
}
dev.off()

pdf("HR_ResidualAutocorrelation.pdf", width = 12, height = 8)
for(i in 1:207){
  if(sum(is.na(missing_normalized_3@data[,3,i])) == 0){
    my_residuals <- missing_normalized_3@data[,3,i] - res$est[,3,i]
    
    acf(my_residuals, 
        na.action = na.pass, main = paste("Subject ", i, ": HR Residual Autocorrelation", sep=""),
        lag.max = 24) 
  }
}
dev.off()


# Subject 207 as an example
pdf("ACF_DBP_Subject207.pdf", width = 8, height = 5)
my_residuals <- missing_normalized_3@data[,1,207] - res$est[,1,207]

acf(my_residuals, 
    na.action = na.pass, main = "DBP: Residual Autocorrelation",
    lag.max = 24) 
dev.off()

my_residuals <- missing_normalized_3@data[,2,207] - res$est[,2,207]

acf(my_residuals, 
    na.action = na.pass, main = "SBP: Residual Autocorrelation",
    lag.max = 24) 

my_residuals <- missing_normalized_3@data[,3,207] - res$est[,3,207]

acf(my_residuals, 
    na.action = na.pass, main = "HR: Residual Autocorrelation",
    lag.max = 24) 


