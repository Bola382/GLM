# ===================================================================================
# Ajustando um MLG com resposta binomial negativa
# ===================================================================================
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())

dados = foreign::read.dta("https://stats.idre.ucla.edu/stat/stata/dae/nb_data.dta")

head(dados)

dados = dados[,-1]

str(dados)

dados$prog = factor(dados$prog)

summary(dados)

# -------------------------------------------------------------------------
# Modelo
# -------------------------------------------------------------------------

modelo = MASS::glm.nb(daysabs ~ math + gender*prog, data = dados, maxit = 300)
design.X = model.matrix(modelo)
coef.fit = coef(modelo)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Redução de covariaveis
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# selecao de variaveis
stepAIC(modelo, direction = "both")

modelo = MASS::glm.nb(daysabs ~ math + gender + prog, data = dados)
design.X = model.matrix(modelo)
coef.fit = coef(modelo)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Diagnostico
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n = nrow(dados) # amostras
p = length(coef.fit) # numero de parametros

# estimativa da media g^(-1)(eta.hat)
mu.hat = exp(design.X%*%coef.fit)

# Estatistica de Pearson
sum(residuals(modelo, type = "pearson")^2)

qchisq(0.05, df = n-p, lower.tail = F) # ok

# Deviance
deviance(modelo)

qchisq(0.05, df = n-p, lower.tail = F) # ok

# **************************
# Residuos
# **************************

library(emulator)
W.matrix = diag(c(modelo$theta*mu.hat/(modelo$theta+mu.hat)))

aux = crossprod(design.X,sqrt(W.matrix))
aux2 = quad.form(W.matrix,design.X)

hat.matrix = quad.form.inv(aux2,aux)

# de Pearson
pearsonRes = residuals(modelo, type="pearson")/sqrt(1-diag(hat.matrix))

# Deviance
devianceRes = residuals(modelo, type = "deviance")/sqrt(1-diag(hat.matrix))

# Quantilicos
qRes = sapply(1:(nrow(dados)), function(a) qnorm(pnbinom(dados$daysabs[a], mu = fitted(modelo)[a], size = modelo$theta)))

# *******************
# envelope
# *******************
B = 1000
qresBoot = matrix(NA, nrow = B, ncol = nrow(dados))
pearsonResBoot = devianceResBoot = qresBoot 

set.seed(2)
for(i in 1:B){
 repl = MASS::rnegbin(n, fitted(modelo),modelo$theta)
 dddd = cbind.data.frame(repl,design.X[,-1])
 mmmm = MASS::glm.nb(repl ~ ., data = dddd, maxit=300)
 
 W.matrix = diag(mmmm$theta*fitted(mmmm)/(mmmm$theta+fitted(mmmm)))
 
 aux = crossprod(design.X,sqrt(W.matrix))
 aux2 = quad.form(W.matrix,design.X)
 
 hat.matrixBoot = quad.form.inv(aux2,aux)
 
 qresBoot[i,] = sort(sapply(1:(nrow(dados)), function(a) qnorm(pnbinom(dddd$repl[a], size = mmmm$theta, mu = fitted(mmmm)[a]))))
 pearsonResBoot[i,] = sort(residuals(mmmm,type = "pearson")/sqrt(1-diag(hat.matrixBoot)))
 devianceResBoot[i,] = sort(residuals(mmmm,type = "deviance")/sqrt(1-diag(hat.matrixBoot)))
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
plot(dados$gender,qRes, pch = 16, cex = .6, 
     ylab = "Resíduo quantílico", xlab = expression(hat(eta)))

# nota de matematica
plot(dados$math,qRes, pch = 16, cex = .6, 
     ylab = "Resíduo quantílico", xlab = expression(hat(eta)));lines(smooth.spline(dados$math,qRes))

# programa
plot(dados$prog,qRes, pch = 16, cex = .6, 
     ylab = "Resíduo quantílico", xlab = expression(hat(eta)))

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

# modelo poisson tbm eh horrivel