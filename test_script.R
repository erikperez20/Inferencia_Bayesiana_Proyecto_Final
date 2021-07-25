# install.packages("interp")

library(fields)
# library(interp)
library(akima)
library(mvtnorm)
filename <- "C:/Users/erik_/Documents/erik documents/Programas/Inferencia_Bayesiana_Proyecto_Final/INFERENCIA_TrueStandardOscillation.txt"

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
# interp(x = delta_array2,y = ldm_array2, z = chi2_array2, xo = 54.7, yo = 2.4356, output = "points")
z0 = interp(x =  delta_array2, y = ldm_array2 , z=chi2_array2, xo = 54.7, yo = 2.4356, extrap = TRUE)
z0

#############################

# Parámetros para las distribuciones a priori:
# Delta cp:
mu_delta <- -165
sd_deltacp <- 51*3 # 3 sigmas de error
# ldm 10^-3:
mu_ldm <- 2.514
sd_ldm <- 0.028*3 # 3 sigmas de error

# Limites del espacio parametrico para cada parametro

# Delta cp:
min_delta <- min(delta_array)
max_delta <- max(delta_array)

min_ldm <- min(ldm_array)
max_ldm <- max(ldm_array)

##############################

posteriori_function <- function(theta){
  # Funcion para calcular la posteriori:
  # theta: vector de parametros
  
  # Distribucion a priori: normal bivariadad con medias 0 y sigmas 1
  priori <- dmvnorm(theta, mean = c(mu_delta,mu_ldm), sigma = matrix(c(sd_deltacp^2,0,0,sd_ldm^2),2,2))
  
  # Likelihood: Exponencial de =chi_cuadrado/2, se utiliza un interpolador debido 
  # a que no se tiene la distribucion conocida de los datos
  z0 <- interp(x =  delta_array2, y = ldm_array2 , z=chi2_array2, xo = theta[1], yo = theta[2])
  likelihood <- exp(-z0$z[1]/2)
  
  return(priori*likelihood)}


## Random Walk proposal distribution:
random_walk_proposal <- function(theta,theta_mean,sj2){
  dist_func <- dmvnorm(theta, mean = theta_mean, sigma = sqrt(sj2)*matrix(c(1,0,0,1),2,2))
  return(dist_func)}


## Candidate generated from the random walk proposal dist.
candidate_random_walk <- function(theta, sj2){
  candidate = rnorm(1, mean = theta, sd = sqrt(sj2))
  return(candidate)}


###### Metropolis-Hasting Algorithm: ######

num_iter = 1000 # number of iterations

# initial sampling values
theta_init <- c(50,2.5) # theta1 = deltacp , theta2 = ldm
p <- length(theta_init)

# tuning parameter for a random walk proposed distribution
s2 = 0.1^2 

##########         Iniciar iteraciones        ##########

theta_val <- theta_init
# Arrays to store the sample values
sample.store <- matrix(0, num_iter, p)

# Start iteration
for(iter in 1:num_iter){
  for (j in 1:p){
    can <- theta_val
    can[j] <- candidate_random_walk(theta_val[j], s2)
    # Agregamos una condicional que limita el espacio paramétrico de los parámetros:
    if (all(can <= c(max_delta,max_ldm) & c(min_delta,min_ldm) <= can & 
        theta_val <= c(max_delta,max_ldm) & c(min_delta,min_ldm) <= theta_val)){
      print('Good')
      # Calculamos el término R
      num1 <- posteriori_function(can) # first term pi(th*|y)
      denom1 <- posteriori_function(theta_val) # second term pi(th^(t-1)|y)
      num2 <- random_walk_proposal(theta_val, can, s2) # 3th term pi(th^(t-1)|th*)
      denom2 <- random_walk_proposal(can, theta_val, s2) #4th term pi(th*|th^(t-1))
      
      # Factor de tasa de aceptación
      AcceptRate <- num1*num2/(denom1*denom2)
      # R factor: minimo entre la tasa de aceptacion y 1
      R <- min(1, AcceptRate)
      print(R)
      # Condition of acceptance of the theta candidate
      if (R == 1){
        theta_val <- can}
      else if (R < 1){
        unif <- runif(1)
        if (unif < R){
          theta_val <- can}
      }
      # Store the sampling values:
      sample.store[iter,] <- can
    }
    else {
      can <- theta_val # si el nuevo candidato no cae en el espacio parametrico regresamos al parametro anterior
      }
    }
}





