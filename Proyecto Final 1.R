# install.packages("interp")
# install.packages("LaplacesDemon")

library(fields)
library("interp")
library(mvtnorm)
library(LaplacesDemon)
filename <- "C:/Users/ALICIA/Documents/Inferencia Bayesiana/INFERENCIA_TrueStandardOscillation.txt"

### Leer el archivo
df <- read.table(file = filename, header = TRUE ,sep = ' ',skip = 1 )
head(df)
### Cambio del nombre de headers
names(df) <- c("delta","m31","chi2")
head(df)

# Arreglo de puntos en el grid para cada par�metro
delta_array <- unique(df$delta)
ldm_array <- unique(df$m31)
chi2_data<- array(df$chi2, dim = c(length(ldm_array),length(delta_array)))
# dim(chi2_data)
image.plot(ldm_array,delta_array, chi2_data)

### Arreglo con todos los valores para la interpolaci�n
delta_array2 <- df$delta
ldm_array2 <- df$m31
chi2_array2 <- df$chi2

#############################

### Par�metros (true) para las distribuciones a priori:
## Delta cp:
mu_delta <- -90
sd_deltacp <- 51*3 # 3 sigmas de error
## ldm/10^-3:
mu_ldm <- 2.514
sd_ldm <- 0.028*3 # 3 sigmas de error

### L�mites del espacio param�trico para cada par�metro
min_delta <- min(delta_array)
max_delta <- max(delta_array)

min_ldm <- min(ldm_array)
max_ldm <- max(ldm_array)

##############################

posteriori_function <- function(theta){
  ### Funci�n para calcular la a posteriori:
  ## theta: vector de par�metros
  
  ## Distribuci�n a priori: normal bivariada con medias valores centrales (true)
  ## y desviaci�n est�ndar 3 sigmas de error
  priori <- dmvnorm(theta, mean = c(mu_delta,mu_ldm), sigma = matrix(c(sd_deltacp^2,0,0,sd_ldm^2),2,2))
  
  # Likelihood: Exponencial de -chi_cuadrado/2, se utiliza un interpolador debido 
  # a que no se tiene la distribucion conocida de los datos
  z0 <- interp(x = delta_array2,y = ldm_array2, z = chi2_array2, xo = theta[1], yo = theta[2], input = "points", output = "points")

  likelihood <- exp(-z0$z[1]/2)
  
  return(priori*likelihood)}


## Random Walk - distribuci�n propuesta:
random_walk_proposal <- function(theta,theta_mean,sj2){
  dist_func <- dmvnorm(theta, mean = theta_mean, sigma = (sj2)*matrix(c(1,0,0,1),2,2))
  return(dist_func)}


## Candidato generado del random walk proposal dist.
candidate_random_walk <- function(theta, sj2){
  candidate = rnorm(1, mean = theta, sd = sqrt(sj2))
  return(candidate)}


###### Metropolis-Hasting Algorithm: ######

num_iter = 1000 # N�mero de iteraciones
burn = 50 # Fijamos el n�mero de muestras del burn-in

## Valores iniciales para el muestreo
theta_init <- c(50,2.5) # theta1 = deltacp , theta2 = ldm
p <- length(theta_init)

## Tuning parameter para el random walk proposal distribution
s2 <- c(5^2,0.1^2)

##########         Iniciar iteraciones        ##########

theta_val <- theta_init
## Arreglos para guardar los valores del muestreo
sample.store <- matrix(0, num_iter, p)

## Inicio de iteraciones
for(iter in 1:num_iter){
  for (j in 1:p){
    can <- theta_val
    can[j] <- candidate_random_walk(theta_val[j], s2[j])
    # Agregamos una condicional que limita el espacio param�trico de los par�metros:
    if (all(can <= c(max_delta,max_ldm) & c(min_delta,min_ldm) <= can & 
            theta_val <= c(max_delta,max_ldm) & c(min_delta,min_ldm) <= theta_val)){
     
      # Calculamos el t�rmino R
      num1 <- posteriori_function(can) # first term pi(th*|y)
      denom1 <- posteriori_function(theta_val) # second term pi(th^(t-1)|y)
      num2 <- random_walk_proposal(theta_val, can, s2) # 3th term pi(th^(t-1)|th*)
      denom2 <- random_walk_proposal(can, theta_val, s2) #4th term pi(th*|th^(t-1))
      
      # Factor de tasa de aceptaci�n. Incluimos los proposal distributions en caso
      # se requiera cambiar a distribuciones no sim�tricas
      AcceptRate <- num1*num2/(denom1*denom2)
      # R factor: m�nimo entre la tasa de aceptaci�n y 1
      R <- min(1, AcceptRate)
    
      # Condiciones de aceptaci�n del candidato
      if (R == 1){
        theta_val <- can}
      else if (R < 1){
        unif <- runif(1)
        if (unif < R){
          theta_val <- can}
      }
      # Guardar los valores del muestreo
      sample.store[iter,] <- can
    }
    else {
      can <- theta_val # Si el nuevo candidato no cae en el espacio param�trico,
      #regresamos al par�metro anterior
    }
  }
  print(iter)
  # Gr�ficos
  if(iter%%10000==0){
    plot(sample.store[1:iter,1],pch=".",main="mu")
    lines(sample.store[1:iter,1])
    plot(sample.store[1:iter,2],pch=".",main="mu")
    lines(sample.store[1:iter,2])
  }
}
### Gr�ficos de la distribuci�n de cada par�metro
d_delta <- density(sample.store[burn:iter,1])
plot(d_delta)
d_ldm <- density(sample.store[burn:iter,2])
plot(d_ldm) 

### Evaluamos la convergencia
Geweke.Diagnostic(sample.store[burn:iter,1])
Geweke.Diagnostic(sample.store[burn:iter,2])

### Para evaluar si se recuperan estimados razonables de los par�metros
sprintf("True mean = %f Delta Posteriori mean = %f",mu_delta,mean(sample.store[burn:iter,1]))
sprintf("True ldm = %f ldm Posteriori mean = %f",mu_ldm,mean(sample.store[burn:iter,2]))

