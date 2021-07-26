install.packages("interp")

library(fields)
library("interp")
filename <- "C:/Users/ALICIA/Documents/Inferencia Bayesiana/INFERENCIA_TrueStandardOscillation.txt"

# Read file
df <- read.table(file = filename, header = TRUE ,sep = ' ',skip = 1 )
head(df)
# Cambio del nombre de headers
names(df) <- c("delta","m31","chi2")
head(df)

# Array of grid points for each parameter
delta_array <- unique(df$delta)
ldm_array <- unique(df$m31)

chi2_data<- array(df$chi2, dim = c(length(ldm_array),length(delta_array)))
dim(chi2_data)

image.plot(ldm_array,delta_array, chi2_data)

delta_array2 <- df$delta
ldm_array2 <- df$m31
chi2_array2 <- df$chi2

# Interpolation: xo and yo are the points where we want to get z
interp(x = delta_array2,y = ldm_array2, z = chi2_array2, xo = 54.76, yo = 2.43567, input = "points", output = "points")
