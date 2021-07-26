library(fields)
library("interp")
library(mvtnorm)
library(LaplacesDemon)


##### DATA
filename <- "C:/Users/ALICIA/Documents/Inferencia Bayesiana/INFERENCIA_TrueStandardOscillation.txt"

### Leer el archivo
df <- read.table(file = filename, header = TRUE ,sep = ' ',skip = 1 )
head(df)
### Cambio del nombre de headers
names(df) <- c("delta","m31","chi2")
head(df)

# Arreglo de puntos en el grid para cada parámetro
delta_array <- unique(df$delta)
ldm_array <- unique(df$m31)
chi2_data<- array(df$chi2, dim = c(length(ldm_array),length(delta_array)))
# dim(chi2_data)
image.plot(ldm_array,delta_array, chi2_data)



### Datos para los ajustes (data generada con otras hipótesis)
filenameH1 <- "C:/Users/ALICIA/Documents/Inferencia Bayesiana/INFERENCIA_TrueStandardOscillation.txt"
dfH1 <- read.table(file = filenameH1, header = TRUE ,sep = ' ',skip = 1 )
head(dfH1)
names(dfH1) <- c("delta","m31","chi2")
head(dfH1)

filenameH2 <- "C:/Users/ALICIA/Documents/Inferencia Bayesiana/INFERENCIA_TrueStandardOscillation.txt"

### Leer el archivo
dfH2 <- read.table(file = filenameH2, header = TRUE ,sep = ' ',skip = 1 )
head(dfH2)
### Cambio del nombre de headers
names(dfH2) <- c("delta","m31","chi2")
head(dfH2)

### Arreglo con todos los valores para la interpolación
delta_array2H1H2 <- dfH1$delta # mismos valores de delta
ldm_array2H1 <- dfH1$m31 # ldm positivo
chi2_array2H1 <- dfH1$chi2

ldm_array2H2 <- dfH2$m31 # ldm negativo
chi2_array2H2 <- dfH2$chi2

#############################

### Parámetros (true) para las distribuciones a priori:
## Delta cp:
mu_delta <- -90
sd_deltacp <- 51*3 # 3 sigmas de error
## ldm/10^-3:
mu_ldm <- 2.514
mu_ldm_inv <- -2.497
sd_ldm <- 0.028*3 # 3 sigmas de error

### Límites del espacio paramétrico para cada parámetro
min_delta <- min(delta_array)
max_delta <- max(delta_array)

min_ldm <- min(ldm_array2H1)
max_ldm <- max(ldm_array2H1)

min_ldminv <- min(ldm_array2H2)
max_ldminv <- max(ldm_array2H2)
##############################

posteriori_function <- function(theta,model){
  ### Función para calcular la a posteriori:
  ## theta: vector de parámetros
  ## model: ldm positivo (1) o negativo (-1)
  
  if(model==1){
  ## Distribución a priori: normal bivariada con medias valores centrales (true)
  ## y desviación estándar 3 sigmas de error
  priori <- dmvnorm(theta, mean = c(mu_delta,mu_ldm), sigma = matrix(c(sd_deltacp^2,0,0,sd_ldm^2),2,2))
  
  # Likelihood: Exponencial de -chi_cuadrado/2, se utiliza un interpolador debido 
  # a que no se tiene la distribucion conocida de los datos
  z0 <- interp(x = delta_array2H1H2,y = ldm_array2H1, z = chi2_array2H1, xo = theta[1], yo = theta[2], input = "points", output = "points")
  }
  else{
    ## Distribución a priori: normal bivariada con medias valores centrales (true)
    ## y desviación estándar 3 sigmas de error
    priori <- dmvnorm(theta, mean = c(mu_delta,mu_ldm_inv), sigma = matrix(c(sd_deltacp^2,0,0,sd_ldm^2),2,2))
    
    # Likelihood: Exponencial de -chi_cuadrado/2, se utiliza un interpolador debido 
    # a que no se tiene la distribucion conocida de los datos
    z0 <- interp(x = delta_array2H1H2,y = ldm_array2H2, z = chi2_array2H2, xo = theta[1], yo = theta[2], input = "points", output = "points")
  }
  likelihood <- exp(-z0$z[1]/2)
  
  return(priori*likelihood)}


## Random Walk - distribución propuesta:
random_walk_proposal <- function(theta,theta_mean,sj2){
  dist_func <- dmvnorm(theta, mean = theta_mean, sigma = (sj2)*matrix(c(1,0,0,1),2,2))
  return(dist_func)}


## Candidato generado del random walk proposal dist.
candidate_random_walk <- function(theta, sj2){
  candidate = rnorm(1, mean = theta, sd = sqrt(sj2))
  return(candidate)}


###### Metropolis-Hasting Algorithm: ######

num_iter = 1000 # Número de iteraciones
burn = 50 # Fijamos el número de muestras del burn-in
model= 1 

## Valores iniciales para el muestreo
theta_init <- c(50,2.5) # theta1 = deltacp , theta2 = ldm
p <- length(theta_init)

## Tuning parameter para el random walk proposal distribution
s2 <- c(5^2,0.1^2) 

##########         Iniciar iteraciones  - ldm positivo      ##########

theta_val <- theta_init
## Arreglos para guardar los valores del muestreo
sample.store <- matrix(0, num_iter, p)

## Inicio de iteraciones
for(iter in 1:num_iter){
  for (j in 1:p){
    can <- theta_val
    can[j] <- candidate_random_walk(theta_val[j], s2[j])
    # Agregamos una condicional que limita el espacio paramétrico de los parámetros:
    if (all(can <= c(max_delta,max_ldm) & c(min_delta,min_ldm) <= can & 
            theta_val <= c(max_delta,max_ldm) & c(min_delta,min_ldm) <= theta_val)){
      
      # Calculamos el término R
      num1 <- posteriori_function(can,model) # first term pi(th*|y)
      denom1 <- posteriori_function(theta_val,model) # second term pi(th^(t-1)|y)
      num2 <- random_walk_proposal(theta_val, can, s2) # 3th term pi(th^(t-1)|th*)
      denom2 <- random_walk_proposal(can, theta_val, s2) #4th term pi(th*|th^(t-1))
      
      # Factor de tasa de aceptación. Incluimos los proposal distributions en caso
      # se requiera cambiar a distribuciones no simétricas
      AcceptRate <- num1*num2/(denom1*denom2)
      # R factor: mínimo entre la tasa de aceptación y 1
      R <- min(1, AcceptRate)
      
      # Condiciones de aceptación del candidato
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
      can <- theta_val # Si el nuevo candidato no cae en el espacio paramétrico,
      #regresamos al parámetro anterior
    }
  }
}

# Gráficos
plot(sample.store[1:iter,1],pch=".",main="delta CP")
lines(sample.store[1:iter,1])
plot(sample.store[1:iter,2],pch=".",main="ldm")
lines(sample.store[1:iter,2])

### Gráficos de la distribución de cada parámetro
d_delta <- density(sample.store[burn:iter,1])
plot(d_delta)
d_ldm <- density(sample.store[burn:iter,2])
plot(d_ldm) 

### Evaluamos la convergencia
Geweke.Diagnostic(sample.store[burn:iter,1])
Geweke.Diagnostic(sample.store[burn:iter,2])

##########         Iniciar iteraciones  - ldm negativo     ##########
model = -1
theta_init <- c(50,-2.5)
theta_val <- theta_init
## Arreglos para guardar los valores del muestreo
sample.store2 <- matrix(0, num_iter, p)

## Inicio de iteraciones
for(iter in 1:num_iter){
  for (j in 1:p){
    can <- theta_val
    can[j] <- candidate_random_walk(theta_val[j], s2[j])
    # Agregamos una condicional que limita el espacio paramétrico de los parámetros:
    if (all(can <= c(max_delta,max_ldminv) & c(min_delta,min_ldminv) <= can & 
            theta_val <= c(max_delta,max_ldminv) & c(min_delta,min_ldminv) <= theta_val)){
      
      # Calculamos el término R
      num1 <- posteriori_function(can,model) # first term pi(th*|y)
      denom1 <- posteriori_function(theta_val,model) # second term pi(th^(t-1)|y)
      num2 <- random_walk_proposal(theta_val, can, s2) # 3th term pi(th^(t-1)|th*)
      denom2 <- random_walk_proposal(can, theta_val, s2) #4th term pi(th*|th^(t-1))
      
      # Factor de tasa de aceptación. Incluimos los proposal distributions en caso
      # se requiera cambiar a distribuciones no simétricas
      AcceptRate <- num1*num2/(denom1*denom2)
      # R factor: mínimo entre la tasa de aceptación y 1
      R <- min(1, AcceptRate)
      
      # Condiciones de aceptación del candidato
      if (R == 1){
        theta_val <- can}
      else if (R < 1){
        unif <- runif(1)
        if (unif < R){
          theta_val <- can}
      }
      # Guardar los valores del muestreo
      sample.store2[iter,] <- can
    }
    else {
      can <- theta_val # Si el nuevo candidato no cae en el espacio paramétrico,
      #regresamos al parámetro anterior
    }
  }

}
# Gráficos
  plot(sample.store2[1:iter,1],pch=".",main="delta CP")
  lines(sample.store2[1:iter,1])
  plot(sample.store2[1:iter,2],pch=".",main="ldm")
  lines(sample.store2[1:iter,2])
  
### Gráficos de la distribución de cada parámetro
d_delta2 <- density(sample.store2[burn:iter,1])
plot(d_delta2)
d_ldm2 <- density(sample.store2[burn:iter,2])
plot(d_ldm2) 

### Evaluamos la convergencia
Geweke.Diagnostic(sample.store2[burn:iter,1])
Geweke.Diagnostic(sample.store2[burn:iter,2])

### Selección de modelos
WAIC1_1 = log(1/mean(sample.store[burn:iter,1])^2)+2*var(log(sample.store[burn:iter,1]))
WAIC1_2 = log(1/mean(sample.store[burn:iter,2])^2)+2*var(log(sample.store[burn:iter,2]))

WAIC2_1 = log(1/mean(sample.store2[burn:iter,1])^2)+2*var(log(sample.store2[burn:iter,1]))
WAIC2_2 = log(1/mean(sample.store2[burn:iter,2])^2)+2*var(log(sample.store2[burn:iter,2]))
