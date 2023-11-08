# ===================================================================================
# Ajustando um MLG com resposta exponencial
# ===================================================================================
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())

data(reliability, package="survival")
dados = ifluid

head(dados)

str(dados)

dados$voltage = factor(dados$voltage)

# breve analise exploratoria

summary(dados)

# -------------------------------------------------------------------------
# Modelo
# -------------------------------------------------------------------------
source("Funcoes auxiliares/escore de Fisher MLG exponencial.R")
aux = lm(time~voltage,data=dados)
design.X = model.matrix(aux)
rm(aux)
coef.fit = FSexp(dados$time,design.X)

coef.fit

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Diagnostico
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n = nrow(dados) # amostras
p = length(coef.fit) # numero de parametros

# estimativa da media g^(-1)(eta.hat)
mu.hat = exp(design.X%*%coef.fit)

# Estatistica de Pearson
sum((dados[,1]-mu.hat)^2/mu.hat^2)

qchisq(0.05, df = n-p, lower.tail = F) # ok

sum((dados[,1]-mu.hat)^2/mu.hat^2)/(n-p) # estimativa de a(phi)

# Deviance
source("Funcoes auxiliares/desvio exponencial.R")
expDev(dados[,1],mu.hat)

qchisq(0.05, df = n-p, lower.tail = F) # ok

expDev(dados[,1],mu.hat)/(n-p) # estimativa de a(phi)

# **************************
# Residuos
# **************************
library(emulator)
aux = crossprod(design.X)

hat.matrix = quad.tform.inv(aux,design.X)

# de Pearson
source("Funcoes auxiliares/residuos de Pearson exponencial.R")

pearsonRes = pearsonResiduals(dados[,1],mu.hat)/sqrt(1-diag(hat.matrix))

# Deviance
devianceRes = expDev(dados[,1],mu.hat,comp=T)/sqrt(1-diag(hat.matrix))

# Quantilicos
rate.hat = 1/mu.hat
qRes = sapply(1:(nrow(dados)), function(a) qnorm(pexp(dados[a,1], rate.hat[a])))

# *******************
# envelope
# *******************
source("Funcoes auxiliares/Simulacao exponencial.R")

B = 1000
qresBoot = matrix(NA, nrow = B, ncol = n)
pearsonResBoot = devianceResBoot = qresBoot 

set.seed(5)
for(i in 1:B){
 repl = rregexp(design.X,coef.fit)
 coef.Boot = FSexp(repl,design.X,print=F)
 
 aux = crossprod(design.X)
 hat.matrixBoot = quad.tform.inv(aux,design.X)
 
 mu.Boot = exp(design.X%*%coef.Boot)
 qresBoot[i,] = sort(sapply(1:n, function(a) qnorm(pexp(repl[a],1/mu.Boot[a]))))
 pearsonResBoot[i,] = sort(pearsonResiduals(repl,mu.Boot)/sqrt(1-diag(hat.matrixBoot)))
 devianceResBoot[i,] = sort(expDev(repl,mu.Boot,comp = T)/sqrt(1-diag(hat.matrixBoot)))
};beepr::beep()
# QQplot residuos de Pearson
bands = apply(pearsonResBoot,2, quantile, c(0.025,.975))
bands2 = apply(pearsonResBoot,2, range)

pearsonRes.sorted = sort(pearsonRes)

color = sapply(1:n, function(a) ifelse(pearsonRes.sorted[a] > bands2[2,a] | pearsonRes.sorted[a] < bands2[1,a], "red", ifelse(pearsonRes.sorted[a] > bands[2,a] | pearsonRes.sorted[a] < bands[1,a], "darkorange", "black")))
size = sapply(color, function(a) ifelse(a == "black", .6, .8))

plot(qnorm((1:n-.5)/n),pearsonRes.sorted,
     xlab = "Quantis N(0,1)", ylab = "Resíduo de Pearson",
     ylim = c(min(bands2[1,]),max(bands[2,])),col=color,
     pch = 16, cex=size)
lines(qnorm((1:n-.5)/n) ,bands[1,], col = "gray50")
lines(qnorm((1:n-.5)/n) ,bands[2,], col = "gray50")
lines(qnorm((1:n-.5)/n), colMeans(pearsonResBoot), lty = 2)
lines(qnorm((1:n-.5)/n) ,bands2[1,])
lines(qnorm((1:(nrow(dados))-.5)/nrow(dados)) ,bands2[2,])
legend("topleft", legend = c("Range", "IC 95%", "Média"), title = "Resíduo simulado", lty = c(1,1,2), col = c("black", "gray50", "black"))


# QQplot residuos de componentes do desvio
bands = apply(devianceResBoot,2, quantile, c(0.025,.975))
bands2 = apply(devianceResBoot,2, range)

devianceRes.sorted = sort(devianceRes)

color = sapply(1:n, function(a) ifelse(devianceRes.sorted[a] > bands2[2,a] | devianceRes.sorted[a] < bands2[1,a], "red", ifelse(devianceRes.sorted[a] > bands[2,a] | devianceRes.sorted[a] < bands[1,a], "darkorange", "black")))
size = sapply(color, function(a) ifelse(a == "black", .6, .8))

plot(qnorm((1:n-.5)/n),devianceRes.sorted,
     xlab = "Quantis N(0,1)", ylab = "Componentes do desvio",
     ylim = c(min(bands2[1,]),max(bands[2,])),col=color,
     pch = 16, cex=size)
lines(qnorm((1:n-.5)/n) ,bands[1,], col = "gray50")
lines(qnorm((1:n-.5)/n) ,bands[2,], col = "gray50")
lines(qnorm((1:n-.5)/n), colMeans(devianceResBoot), lty = 2)
lines(qnorm((1:n-.5)/n) ,bands2[1,])
lines(qnorm((1:(nrow(dados))-.5)/nrow(dados)) ,bands2[2,])
legend("topleft", legend = c("Range", "IC 95%", "Média"), title = "Resíduo simulado", lty = c(1,1,2), col = c("black", "gray50", "black"))

# QQplot residuos quantilicos
bands = apply(qresBoot,2, quantile, c(0.025,.975))
bands2 = apply(qresBoot,2, range)

qRes.sorted = sort(qRes)

color = sapply(1:n, function(a) ifelse(qRes.sorted[a] > bands2[2,a] | qRes.sorted[a] < bands2[1,a], "red", ifelse(qRes.sorted[a] > bands[2,a] | qRes.sorted[a] < bands[1,a], "darkorange", "black")))
size = sapply(color, function(a) ifelse(a == "black", .6, .8))

plot(qnorm((1:n-.5)/n),qRes.sorted,
     xlab = "Quantis N(0,1)", ylab = "Resíduo quantilico",
     ylim = c(min(bands2[1,]),max(bands[2,])),
     col=color,
     pch = 16, cex=size)
lines(qnorm((1:n-.5)/n) ,bands[1,], col = "gray50")
lines(qnorm((1:n-.5)/n) ,bands[2,], col = "gray50")
lines(qnorm((1:n-.5)/n), colMeans(qresBoot), lty = 2)
lines(qnorm((1:n-.5)/n) ,bands2[1,])
lines(qnorm((1:(nrow(dados))-.5)/nrow(dados)) ,bands2[2,])
legend("topleft", legend = c("Range", "IC 95%", "Média"), title = "Resíduo simulado", lty = c(1,1,2), col = c("black", "gray50", "black"))

# ***********************
# residuos X mu.hat
# ***********************

# residuo de Pearson
plot(mu.hat,pearsonRes, pch = 16, cex = .6, xlab = expression(hat(mu)), ylab = "Residuo de Pearson")
lines(smooth.spline(mu.hat,pearsonRes))

# componentes do desvio
plot(mu.hat,devianceRes, pch = 16, cex = .6, xlab = expression(hat(mu)), ylab = "Componentes do desvio")
lines(smooth.spline(mu.hat,devianceRes))

# residuo quantilico
plot(mu.hat,qRes, pch = 16, cex = .6, xlab = expression(hat(mu)), ylab = "Resíduo quantilico")
lines(smooth.spline(mu.hat,qRes))

# ***************************
# Residuos x preditor linear
# ***************************

# Pearson
plot(design.X%*%coef.fit,pearsonRes, pch = 16, cex = .6, 
     ylab = "Resíduos de Pearson", xlab = expression(hat(eta)));lines(smooth.spline(design.X%*%coef.fit,pearsonRes))

# Componentes do desvio
plot(design.X%*%coef.fit,devianceRes, pch = 16, cex = .6, 
     ylab = "Componentes do desvio", xlab = expression(hat(eta)));lines(smooth.spline(design.X%*%coef.fit,devianceRes))

# Quantilico
plot(design.X%*%coef.fit,qRes, pch = 16, cex = .6, 
     ylab = "Resíduo quantílico", xlab = expression(hat(eta)));lines(smooth.spline(design.X%*%coef.fit,qRes))

# ***************************
# Residuos x covariavel
# ***************************

plot(qRes~dados$voltage, pch = 16, cex = .6, 
     ylab = "Resíduo quantílico", xlab = expression(hat(eta)))

# ---------------------------------------------------------------------------
# Pontos influentes
# ---------------------------------------------------------------------------

plot(diag(hat.matrix),pch=16,cex=.6,xlab="",ylab="Alavancagem",
     col = ifelse(diag(hat.matrix)>2*p/n,"red","black"));abline(h=2*p/n)
