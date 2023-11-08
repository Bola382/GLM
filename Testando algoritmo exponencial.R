# ===================================================================================
# Verificando a implementacao do algoritmo Escore de Fisher com dados simulados 
# ===================================================================================
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())

source("funcoes auxiliares/Simulacao exponencial.R")
source("funcoes auxiliares/escore de Fisher MLG exponencial.R")
source("funcoes auxiliares/desvio exponencial.R")
# ----------------------
# Gerando dados
# ----------------------

set.seed(1)
n = 500
X = cbind(1, rnorm(n, mean = 0, sd = 1), rgamma(n, shape=5, rate=2))
beta = c(2, 2, -.5)

y = rregexp(X, beta)

plot(X[,2],y, pch =16, cex=.6, xlab = expression(x[1]))
plot(X[,3],y, pch =16, cex=.6, xlab = expression(x[2]))

coef.fit = FSexp(y,X, print = T)

coef.fit

dados = cbind(y,X[,-1])

# inversa da informacao de Fisher esperada
K.inv = solve(crossprod(X))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MCMC
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
R = 1000
set.seed(1)
MC = replicate(R,rregexp(X, beta))
betaMC = apply(MC,2,function(a) FSexp(a,X,print=F))

# teoricamente beta.hat ~ N(beta, K.inv)
hist(betaMC[1,], freq=F,breaks=16, ylim = c(0,4.5), 
     main = "", ylab = "Densidade",xlab=expression(hat(beta)[0]))  
curve(dnorm(x,mean=beta[1],sd = sqrt(K.inv[1,1])),add=T)

hist(betaMC[2,], freq=F,breaks=16, ylim = c(0,10), 
     main = "", ylab = "Densidade",xlab=expression(hat(beta)[1]))  
curve(dnorm(x,mean=beta[2],sd = sqrt(K.inv[2,2])),add=T)

hist(betaMC[3,], freq=F, ylim = c(0,10), 
     main = "", ylab = "Densidade",xlab=expression(hat(beta)[2]))  
curve(dnorm(x,mean=beta[3],sd = sqrt(K.inv[3,3])),add=T)

#---------------------------------------------
# verificando taxa de cobertura dos intervalos
#---------------------------------------------

LI = betaMC[3,]-1.96*sqrt(K.inv[3,3])
LS = betaMC[3,]+1.96*sqrt(K.inv[3,3])
plot(1, xlim = c(1,R), ylim = c(min(LI)-.1,max(LS)+.1), 
     xlab="Replicação", ylab = "",type ="n")
c=0
for(i in 1:R){
 if(beta[3] > LI[i] & beta[3] < LS[i]){
  c = c + 1;color="gray"
  }else{color="black"}
 segments(x0 = i, x1 = i, y0 = LI[i], y1 = LS[i], 
          col = color)
}
abline(h=beta[3], col=2, lty=2)
text(150,-.2,labels=paste0("Taxa de cobertura: ", round(c/R*100,2),"%"))
legend("bottomleft",legend = expression(beta[2]), lty=2,col=2, bty = "n")

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
}
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

# ***********************
# Correlacao serial
# ***********************

# Pearson
color = sapply(1:n, function(a) ifelse(pearsonRes[a] > 3 | pearsonRes[a] < -3, "red", ifelse(pearsonRes[a] > 2 | pearsonRes[a] < -2, "darkorange", "black")))
size = sapply(color, function(a) ifelse(a == "black", .6, .8))

plot(pearsonRes, ylim = c(-3,3)+c(-.1,+.5), pch = 16, cex = size, 
     ylab = "Resíduos de Pearson", xlab = "Indice",
     col = color)
abline(h=mean(qRes),lty=2)
abline(h=c(-3,3), col = "black")
abline(h=c(-2,2), col = "gray50")

# Componentes do desvio
color = sapply(1:n, function(a) ifelse(devianceRes[a] > 3 | devianceRes[a] < -3, "red", ifelse(devianceRes[a] > 2 | devianceRes[a] < -2, "darkorange", "black")))
size = sapply(color, function(a) ifelse(a == "black", .6, .8))

plot(devianceRes, ylim = c(-3,3)+c(-.1,+.1), pch = 16, cex = size, 
     ylab = "Componentes do desvio", xlab = "Indice",
     col = color)
abline(h=mean(qRes),lty=2)
abline(h=c(-3,3), col = "black")
abline(h=c(-2,2), col = "gray50")

# quantilicos
color = sapply(1:n, function(a) ifelse(qRes[a] > 3 | qRes[a] < -3, "red", ifelse(qRes[a] > 2 | qRes[a] < -2, "darkorange", "black")))
size = sapply(color, function(a) ifelse(a == "black", .6, .8))

plot(qRes, ylim = c(-3,3)+c(-.1,+.1), pch = 16, cex = size, 
     ylab = "Resíduos quantilicos", xlab = "Indice",
     col = color)
abline(h=mean(qRes),lty=2)
abline(h=c(-3,3), col = "black")
abline(h=c(-2,2), col = "gray50")

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

# genero
plot(dados[,3],qRes, pch = 16, cex = .6, 
     ylab = "Resíduo quantílico", xlab = expression(hat(eta)));lines(smooth.spline(dados[,3],qRes))

# nota de matematica
plot(dados[,4],qRes, pch = 16, cex = .6, 
     ylab = "Resíduo quantílico", xlab = expression(hat(eta)));lines(smooth.spline(dados[,4],qRes))

# ---------------------------------------------------------------------------
# Pontos influentes
# ---------------------------------------------------------------------------

plot(diag(hat.matrix),pch=16,cex=.6,xlab="",ylab="Alavancagem",
     col = ifelse(diag(hat.matrix)>2*p/n,"red","black"));abline(h=2*p/n)

plot(dffits(modelo),pch=16,cex=.6,col=ifelse(dffits(modelo)>2*sqrt(p/n),"red","black"),
     ylab = "DFFitts", xlab = "")
abline(h=2*sqrt(p/n))

cook.modelo = cooks.distance(modelo)
cook.f = pf(cook.modelo, df1 = p, df2 = n-p, lower.tail = F)
plot(cook.modelo, pch = 16, cex=.6, col = ifelse(cook.f <= .5, "red", "black"))
abline(h = qf(.5,df1=p,df2=n-p))
