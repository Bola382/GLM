# ===================================================================================
# Verificando a implementacao do algoritmo Escore de Fisher com dados simulados 
# ===================================================================================
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())

source("funcoes auxiliares/Simulacao geometrica.R")
source("funcoes auxiliares/escore de Fisher MLG geometrica.R")
source("funcoes auxiliares/desvio geometrica.R")
# ----------------------
# Gerando dados
# ----------------------

set.seed(1)
n = 50
X = cbind(1, rnorm(n, mean = 0, sd = 1), rgamma(n, shape=5, rate=2))
beta = c(1, -2, .5)

y = rreggeom(X, beta)

plot(X[,2],y, pch =16, cex=.6, xlab = expression(x[1]))
plot(X[,3],y, pch =16, cex=.6, xlab = expression(x[2]))

bb = FSgeo(y,X, print = T)

bb$result
bb$historico

aa = glm(y ~ X[,2]+X[,3], family = MASS::negative.binomial(theta=1)) 

bb$result
aa$coef

# ------------------
# Monte Carlo
# ------------------
R = 1000
MC = replicate(R, rreggeom(X, beta))

MCbeta = apply(MC, 2, function(a) FSgeo(a,X,print=F)$result)
MCbeta2 = apply(MC, 2, function(a) glm(a~X[,2]+X[,3],  maxit = 300,family = MASS::negative.binomial(theta = 1))$coef)

# beta 0
par(mar = c(5.1, 4.5, 3, 1))
plot(MCbeta[1,], MCbeta2[1,], cex = .6,pch = 16, 
     xlab = expression(beta[0]^"FS"), ylab = expression(beta[0]^"R"))
abline(a=0,b=1, col = 2, lty = 2, lwd = 2)
legend("topleft", legend = "y = x", lty = 2, col = 2, lwd = 2)       

# beta 1
plot(MCbeta[2,], MCbeta2[2,], cex = .6,pch = 16, 
     xlab = expression(beta[1]^"FS"), ylab = expression(beta[1]^"R"))
abline(a=0,b=1, col = 2, lty = 2, lwd = 2)
legend("topleft", legend = "y = x", lty = 2, col = 2, lwd = 2)       

# beta 2
plot(MCbeta[3,], MCbeta2[3,], cex = .6,pch = 16, 
     xlab = expression(beta[2]^"FS"), ylab = expression(beta[2]^"R"))
abline(a=0,b=1, col = 2, lty = 2, lwd = 2)
legend("topleft", legend = "y = x", lty = 2, col = 2, lwd = 2)       
