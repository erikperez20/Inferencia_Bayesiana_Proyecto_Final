filename <- "C:/Users/erik_/Documents/erik documents/Programas/Inferencia_Bayesiana_Proyecto_Final/INFERENCIA.txt"

# Read file
df <- read.table(file = filename, header = TRUE ,sep = ' ',skip = 1 )
head(df)

# Cambio del nombre de headers
names(df) <- c("delta","theta","m31","chi2")
head(df)

# Array of grid points for each parameter
delta_array <- unique(df$delta)
theta_array <- unique(df$theta)
ldm_array <- unique(df$m31)

chi2_data <- array(df$chi2, dim = c(length(ldm_array),length(theta_array),length(delta_array)))
dim(chi2_data[,,1])

library(fields)
image.plot(ldm_array,theta_array, chi2_data[,,100])
image.plot(ldm_array,delta_array, chi2_data[,24,])
image.plot(theta_array,delta_array, chi2_data[7,,])



# chi2_data[1,,]
# 
# chi2_data[,,1]
# 
# dim(chi2_data[201,,])
# # lapply(list(df$delta), function(x) rle(x)$lengths)
