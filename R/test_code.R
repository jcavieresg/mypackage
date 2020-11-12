rm(list = ls())
setwd("C:/Users/Usuario/Desktop/Thesis/model/experimental_codes")

library("R.matlab")
library(Matrix)

#dsites = R.matlab::readMat("Data2D_25h.mat")
dsites = R.matlab::readMat("Data2D_1089h.mat")

#dsites = R.matlab::readMat("Data2D_25h.mat")
ctrs = R.matlab::readMat("Data2D_81u.mat")
ctrs$dsites

class(dsites) # lista
class(ctrs)   # lista

# Lo convertimos a matriz
dsites = dsites$dsites
ctrs = ctrs$dsites

class(dsites)
class(ctrs)


# Comparamos
par(mfrow=c(1,2))
plot(dsites[, 1], dsites[,2])
plot(ctrs[, 1], ctrs[,2])


library(Rcpp)
#Rcpp::sourceCpp("smooth_TPS.cpp")
Rcpp::sourceCpp("RBF_alterna.cpp")
#Rcpp::sourceCpp("example_RBF_3.cpp")


#===================================================
#         Calculo de matriz de distancias
#
# Obs: dsites = Sitios con observaciones
#      ctrs = Centros para evaluar las distancias
#
# Nota: Generalmente se asume que dsites = ctrs
#===================================================

DM_data = DistanceMatrix(dsites, ctrs)
DM_data2 = DistanceMatrix2(dsites, ctrs)

class(DM_data)
class(DM_data2)

DM_data
DM_data2

#par(mfrow=c(2,2))
heatmap(as.matrix(DM_data))
heatmap(as.matrix(DM_data2))

library(lattice)
levelplot(DM_data, col.regions=heat.colors(100))

library(viridisLite)
coul <- viridis(100)
levelplot(DM_data, col.regions = coul) 


#write.csv(DM_data,'DM_data.csv')

#=====================================================================
#                           Calculos de phi
#=====================================================================
# Funcion para calcular los phi asociados a la matriz
# de distancia calculada

# funcion: radialFunction(1, 2, 3)
# Argumentos: 1) Matriz de distancias (o r de acuerdo a la literatura)
#             2) Tipo de RBF
# Globales
# 1 = Spline lineal
# 2 = Thin Plate Spline
# 3 = Gaussian
# 4 = Multicuadratic
# 5 = Inversa Multicuadratic

# Soporte compacto
# 6 = Wenland 
# 
#
#             3) Valor entero de R para funciones globales o compactas              
#
#_____________________________________________________________________
#                            EJEMPLOS
#                     Definimos los parámetros
#=====================================================================
Lineal = 1
TPS = 2
Gaussian = 3
Multi = 4
R = 1

phi_lineal = radialFunction(DM_data, Lineal, R)

# Thin plate Spline RBF
phi_TPS = radialFunction(DM_data, TPS, R)

# Gaussian RBF
phi_Gaussian = radialFunction(DM_data, Gaussian, R)

# Multicuadratic RBF
phi_Multi = radialFunction(DM_data, Multi, R)
#=====================================================================




#=======================================================================
#                Ajuste mediante minimos cuadrados
#=======================================================================
# Funcion para calcular RBF mediante minimos cuadrados 
# funcion: RBF_LS(1, 2, 3, 4)
# Argumentos: 1 = dsites: matriz 2X2 con posiciones en el plano en forma
#                aleatoria)
#             2 = ctrs: matriz 2X2 con posiciones centradas en el plano 
#                 en forma uniforme
#             3 = Tipo de RBF a utilizar (1 = lineal, 2 = TPS, 3 = Gaussian)
#                 4 = Multicuadratica                 
#             4 = Entero que diferencia entre funciones globales y de soporte
#                 Compacto
#_____________________________________________________________________
#_____________________________________________________________________
#                                    EJEMPLOS

#======================================================
# Ejemplo 1
# Ajuste mediante minimos cuadrados y una función creada 
fit_LS = RBF_LS(dsites, ctrs, 2, 1, 100) # Usando funcion radialFunction
fit_LS$phi_coefficients
fit_LS$Pf
fit_LS$maxerr
fit_LS$Distances

#======================================================
# Ejemplo 2
# Simulo un 1 con distribucion normal y de largo dsites
y = rnorm(length(dsites[,1]), 0, 1)
hist(y)

# Ajuste con una variable respuesta y
fit2_LS = RBF_LS2(y, dsites, ctrs, 2, 1, 100)
fit2_LS$phi_coefficients
fit2_LS$Pf
fit2_LS$maxerr
fit2_LS$Distances




#=======================================================================
#       Smoothing Thin Plate Spline (minimos cuadrados penalizados)
#=======================================================================
# Funcion para calcular RBF mediante minimos cuadrados penalizados 
# funcion: RBF_LSP(1, 2, 3, 4)
# Argumentos: 1 = dsites: matriz 2X2 con posiciones en el plano en forma
#                aleatoria)
#             2 = ctrs: matriz 2X2 con posiciones centradas en el plano 
#                 en forma uniforme
#             3 = Tipo de RBF a utilizar (1 = lineal, 2 = TPS, 3 = Gaussian)
#                 4 = Multicuadratica                 
#             4 = Entero que diferencia entre funciones globales y de soporte
#                 Compacto
#             5 =  Número de evaluaciones
#             6 = omega (controla el smooth en el ajuste)
#_____________________________________________________________________



#=================================================================================
# Ejemplo 1
# RBF penalized TPS
# Aquí uso dsites = ctrs (por codigo de matlab) y funcion creada dentro del codigo
# para evaluacion
omega1 = 1
omega10 = 10
omega100 = 100
omega1000 = 1000

neval = 100

fit_tps1 = RBF_LSP(dsites, dsites, TPS, R, neval, omega1)
fit_tps2 = RBF_LSP(dsites, dsites, TPS, R, neval, omega10)
fit_tps3 = RBF_LSP(dsites, dsites, TPS, R, neval, omega100)
fit_tps4 = RBF_LSP(dsites, dsites, TPS, R, neval, omega1000)                   


length(fit_tps1$Pf)
length(fit_tps1$epoints[,1])

par(mfrow=c(2,2))
library(plot3D)
scatter3D(fit_tps1$epoints[, 1], fit_tps1$epoints[, 2], fit_tps1$Pf, theta = 160, phi = 10, bty = "g")

scatter3D(fit_tps1$epoints[, 1], fit_tps1$epoints[, 2], fit_tps2$Pf, theta = 160, phi = 10, bty = "g")

scatter3D(fit_tps1$epoints[, 1], fit_tps1$epoints[, 2], fit_tps3$Pf, theta = 160, phi = 10, bty = "g")

scatter3D(fit_tps1$epoints[, 1], fit_tps1$epoints[, 2], fit_tps4$Pf, theta = 160, phi = 10, bty = "g")

error_tps1 = fit_tps1$maxerr
error_tps2 = fit_tps2$maxerr
error_tps3 = fit_tps3$maxerr
error_tps4 = fit_tps4$maxerr


table_error = c(error_tps1, error_tps2, error_tps3, error_tps4)
table_error


# Matlab errores maximos 
w1000 = 7.580288e-002 
w100  = 4.215865e-002 
w10   = 3.349371e-002 
w1    = 1.041350e-001
errores_matlab = c(w1000, w100, w10, w1)
errores_matlab



x = 1:100;
y = x^2;
plot(x, y);
plot(log(x), log(y));






#=================================================================================
# Ejemplo 2
# RBF penalized TPS
#                           DATOS SIMULADOS DE dsites y ctrs

dsites_sim = cbind(x1 = runif(4000), x2 = runif(4000))
ctrs_sim =   cbind(x1 = runif(4000), x2 = runif(4000))

dsites_sim = as.matrix(dsites_sim)
ctrs_sim = as.matrix(ctrs_sim)

plot(dsites_sim)
plot(ctrs_sim)

omega1 = 1
omega10 = 10
omega100 = 100
omega1000 = 1000

neval = 100


fit_tps1 = RBF_LSP(dsites_sim, ctrs_sim, TPS, R, neval, omega1)
fit_tps2 = RBF_LSP(dsites_sim, ctrs_sim, TPS, R, neval, omega10)
fit_tps3 = RBF_LSP(dsites_sim, ctrs_sim, TPS, R, neval, omega100)
fit_tps4 = RBF_LSP(dsites_sim, ctrs_sim, TPS, R, neval, omega1000)                   


length(fit_tps1$Pf)
length(fit_tps1$epoints[,1])

par(mfrow=c(2,2))
library(plot3D)
scatter3D(fit_tps1$epoints[, 1], fit_tps1$epoints[, 2], fit_tps1$Pf, theta = 160, phi = 10, bty = "g")

scatter3D(fit_tps1$epoints[, 1], fit_tps1$epoints[, 2], fit_tps2$Pf, theta = 160, phi = 10, bty = "g")

scatter3D(fit_tps1$epoints[, 1], fit_tps1$epoints[, 2], fit_tps3$Pf, theta = 160, phi = 10, bty = "g")

scatter3D(fit_tps1$epoints[, 1], fit_tps1$epoints[, 2], fit_tps4$Pf, theta = 160, phi = 10, bty = "g")

error_tps1 = fit_tps1$maxerr
error_tps2 = fit_tps2$maxerr
error_tps3 = fit_tps3$maxerr
error_tps4 = fit_tps4$maxerr


table_error = c(error_tps1, error_tps2, error_tps3, error_tps4)
table_error










#=================================================================================
# Ejemplo 3
# RBF penalized TPS
# Simulo una variale y ~ N(0,1) para ajuste de datos "reales"
dsites_sim = cbind(x1 = runif(1000), x2 = runif(1000))
ctrs_sim =   cbind(x1 = runif(1000), x2 = runif(1000))

dsites_sim = as.matrix(dsites_sim)
ctrs_sim = as.matrix(ctrs_sim)

plot(dsites_sim)
plot(ctrs_sim)


y = rnorm(length(dsites_sim[,1]), 0, 1)
length(y)

SPL = 1
TPS = 2
Gaussian = 3
Multi = 4

omega1 = 1
omega10 = 10
omega100 = 100
omega1000 = 1000

neval = 50

# 
fit_tps1 = RBF_LSP2(y, dsites_sim, dsites_sim, TPS, 1, neval, omega1000)
fit_tps2 = RBF_LSP2(y, dsites_sim, dsites_sim, TPS, 1, neval, omega100)
fit_tps3 = RBF_LSP2(y, dsites_sim, dsites_sim, TPS, 1, neval, omega10)
fit_tps4 = RBF_LSP2(y, dsites_sim, dsites_sim, TPS, 1, neval, omega1)


par(mfrow=c(2,2))
library(plot3D)
scatter3D(fit_tps1$epoints[, 1], fit_tps1$epoints[, 2], fit_tps1$Pf, theta = 160, phi = 10, bty = "g")

scatter3D(fit_tps1$epoints[, 1], fit_tps1$epoints[, 2], fit_tps2$Pf, theta = 160, phi = 10, bty = "g")

scatter3D(fit_tps1$epoints[, 1], fit_tps1$epoints[, 2], fit_tps3$Pf, theta = 160, phi = 10, bty = "g")

scatter3D(fit_tps1$epoints[, 1], fit_tps1$epoints[, 2], fit_tps4$Pf, theta = 160, phi = 10, bty = "g")



error_tps1 = fit_tps1$maxerr
error_tps2 = fit_tps2$maxerr
error_tps3 = fit_tps3$maxerr
error_tps4 = fit_tps4$maxerr


table_error = c(error_tps1, error_tps2, error_tps3, error_tps4)
table_error







#=================================================================================
# Ejemplo 4
# RBF penalized TPS con SINGULAR VALEUE DESCOMPOSITION (SVD)

# Simulo una variale y ~ N(0,1) para ajuste de datos "reales"
dsites_sim = cbind(x1 = runif(1000), x2 = runif(1000))
ctrs_sim =   cbind(x1 = runif(1000), x2 = runif(1000))

dsites_sim = as.matrix(dsites_sim)
ctrs_sim = as.matrix(ctrs_sim)

plot(dsites_sim)
plot(ctrs_sim)


y = rnorm(length(dsites_sim[,1]), 0, 1)
length(y)

SPL = 1
TPS = 2
Gaussian = 3
Multi = 4

omega1 = 1
omega10 = 10
omega100 = 100
omega1000 = 1000

neval = 50

# 
fit_tps1 = RBF_LSPSVD(y, dsites_sim, dsites_sim, TPS, 1, neval, omega1000)
fit_tps2 = RBF_LSPSVD(y, dsites_sim, dsites_sim, TPS, 1, neval, omega100)
fit_tps3 = RBF_LSPSVD(y, dsites_sim, dsites_sim, TPS, 1, neval, omega10)
fit_tps4 = RBF_LSPSVD(y, dsites_sim, dsites_sim, TPS, 1, neval, omega1)


par(mfrow=c(2,2))
library(plot3D)
scatter3D(fit_tps1$epoints[, 1], fit_tps1$epoints[, 2], fit_tps1$Pf, theta = 160, phi = 10, bty = "g")

scatter3D(fit_tps1$epoints[, 1], fit_tps1$epoints[, 2], fit_tps2$Pf, theta = 160, phi = 10, bty = "g")

scatter3D(fit_tps1$epoints[, 1], fit_tps1$epoints[, 2], fit_tps3$Pf, theta = 160, phi = 10, bty = "g")

scatter3D(fit_tps1$epoints[, 1], fit_tps1$epoints[, 2], fit_tps4$Pf, theta = 160, phi = 10, bty = "g")



error_tps1 = fit_tps1$maxerr
error_tps2 = fit_tps2$maxerr
error_tps3 = fit_tps3$maxerr
error_tps4 = fit_tps4$maxerr


table_error = c(error_tps1, error_tps2, error_tps3, error_tps4)
table_error






#========================================================================
#                    Ajuste mediante otro tipo de RBF
#========================================================================

fit_lineal = RBF_LSP2(y, dsites, dsites, SPL, 1, neval, omega100)
fit_tps = RBF_LSP2(y, dsites, dsites, TPS, 1, neval, omega100)
fit_gaussian = RBF_LSP2(y, dsites, dsites, Gaussian, 1, neval, omega100)
fit_multi = RBF_LSP2(y, dsites, dsites, Multi, 1, neval, omega100)

dim(fit_lineal$epoints)
length(fit_lineal$Pf)


# Creo Gráficos
library(plot3D)
x11()
par(mfrow=c(1,2))
scatter3D(fit_lineal$epoints[, 1], fit_lineal$epoints[, 2], fit_lineal$Pf, theta = 5, phi = 25,
          bty = "g", colkey = FALSE)


scatter3D(fit_tps$epoints[, 1], fit_tps$epoints[, 2], fit_tps$Pf, theta = 5, phi = 25,
          bty = "g", colkey = FALSE)


scatter3D(fit_gaussian$epoints[, 1], fit_gaussian$epoints[, 2], fit_gaussian$Pf, theta = 5, phi = 25,
          bty = "g", colkey = FALSE)

scatter3D(fit_multi$epoints[, 1], fit_multi$epoints[, 2], fit_multi$Pf, theta = 5, phi = 25,
          bty = "g", colkey = FALSE)


error_lineal = fit_lineal$maxerr
error_tps = fit_tps$maxerr
error_gaussian = fit_gaussian$maxerr
error_multi = fit_multi$maxerr


table_error = c(error_lineal, error_tps, error_gaussian, error_multi)
table_error






#====================================================================================
#                               Ajuste de datos reales
#====================================================================================

setwd("C:/Users/Usuario/Desktop/Projects/INLA_isotope/data")


# load data 
dat = read.csv("data_total.csv", header=T)
dim(dat)
head(dat)
str(dat)

options(scipen=99)
hist(dat$N)

head(dat)

dsites = cbind(dat$Lat, dat$Lon)
y = dat$N


tps = RBF_LSP2(y, dsites, dsites, TPS, 1, neval, omega1000)
gaussian = RBF_LSP2(y, dsites, dsites, Gaussian, 1, neval, omega1000)


tps_error = tps$maxerr
gaussian_error = gaussian$maxerr

table_error2 = c(tps_error, gaussian_error)
table_error2

par(mfrow=c(1,2))
scatter3D(tps$epoints[, 1], tps$epoints[, 2], tps$Pf, theta = 5, phi = 25,
          bty = "g", colkey = FALSE)


scatter3D(gaussian$epoints[, 1], gaussian$epoints[, 2], gaussian$Pf, theta = 5, phi = 25,
          bty = "g", colkey = FALSE)


install.packages("radix")
install.packages("rticles")



