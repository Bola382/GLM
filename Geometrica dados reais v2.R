# ===================================================================================
# Ajustando um MLG com resposta geometrica
# ===================================================================================
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())

dados = foreign::read.dta("https://stats.idre.ucla.edu/stat/stata/dae/nb_data.dta")

head(dados)

dados = dados[,-1]

str(dados)

dados$prog = factor(dados$prog)

# breve analise exploratoria

summary(dados)

plot(dados$daysabs~dados$gender)
plot(dados$daysabs,dados$math)
plot(dados$daysabs~dados$prog)

# -------------------------------------------------------------------------
# Modelo
# -------------------------------------------------------------------------

modelo = glm(daysabs ~ math + gender*prog, data = dados,family = MASS::negative.binomial(theta = 1, link = "log"))
design.X = model.matrix(modelo)
coef.fit = coef(modelo)

source("Funcoes auxiliares/escore de Fisher MLG geometrica.R")
modelo.meu = FSgeo(dados$daysabs, design.X)

cbind(coef.fit,modelo.meu$result, abs(coef.fit-modelo.meu$result))
# IC modelo glm
library(emulator)
W.matrix = diag(c(fitted(modelo)/(1+fitted(modelo))))

K.inv = solve(quad.form(W.matrix, design.X))


LI = coef.fit - 1.96*sqrt(diag(K.inv))
LS = coef.fit + 1.96*sqrt(diag(K.inv))

cbind(LI,LS)

# IC modelo meu

mu.hat = exp(design.X%*%modelo.meu$result)

W.matrix = diag(c(mu.hat/(1+mu.hat)))

K.inv = solve(quad.form(W.matrix, design.X))


LI = coef.fit - 1.96*sqrt(diag(K.inv))
LS = coef.fit + 1.96*sqrt(diag(K.inv))

cbind(LI,LS)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Redução de covariaveis
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# correlacao entre variaveis
boxplot(dados$math~dados$gender, names = c("Feminino","Masculino"),xlab = "Gênero", ylab = "Nota de matemática")
boxplot(dados$math~dados$prog, xlab = "Programa", ylab = "Nota de matemática")
table(dados$gender,dados$prog)

cor.test(table(dados$gender,dados$prog)[1,],table(dados$gender,dados$prog)[2,], method = "spearman", exact = T)

# selecao de variaveis
stepAIC(modelo, direction = "both")

modelo = glm(daysabs ~ math + gender + prog, data = dados,family = MASS::negative.binomial(theta = 1, link = "log"))
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
sum((dados$daysabs-mu.hat)^2/(mu.hat+mu.hat^2))
sum(residuals(modelo, type = "pearson")^2)

qchisq(0.05, df = n-p, lower.tail = F) # ok

sum(residuals(modelo, type = "pearson")^2)/(n-p) # estimativa de a(phi)

# Deviance
source("Funcoes auxiliares/desvio geometrica.R")
geoDev(dados$daysabs,mu.hat)
deviance(modelo)
sum(residuals(modelo, type = "deviance")^2)

qchisq(0.05, df = n-p, lower.tail = F) # ok

deviance(modelo)/(n-p) # estimativa de a(phi)

# **************************
# Residuos
# **************************
library(emulator)
W.matrix = diag(c(mu.hat/(1+mu.hat)))

aux = crossprod(design.X,sqrt(W.matrix))
aux2 = quad.form(W.matrix,design.X)

hat.matrix = quad.form.inv(aux2,aux)

# de Pearson
pearsonRes = residuals(modelo, type="pearson")/sqrt(1-diag(hat.matrix))

# Deviance
devianceRes = residuals(modelo, type = "deviance")/sqrt(1-diag(hat.matrix))

# Quantilicos
p.hat = c(1/(mu.hat+1))
qRes = sapply(1:(nrow(dados)), function(a) qnorm(pgeom(dados$daysabs[a], p.hat[a])))

# *******************
# envelope
# *******************
source("Funcoes auxiliares/Simulacao geometrica.R")

B = 1000
qresBoot = matrix(NA, nrow = B, ncol = n)
pearsonResBoot = devianceResBoot = qresBoot 

set.seed(5)
for(i in 1:B){
 repl = rreggeom(design.X,coef.fit)
 dddd = cbind.data.frame(repl,design.X[,-1])
 mmmm = glm(repl ~ ., data = dddd,maxit = 500,family = MASS::negative.binomial(theta = 1, link = "log"))
 
 W.matrix = diag(fitted(mmmm)/(1+fitted(mmmm)))
 
 aux = crossprod(design.X,sqrt(W.matrix))
 aux2 = quad.form(W.matrix,design.X)
 
 hat.matrixBoot = quad.form.inv(aux2,aux)
 
 mu = c(1/(exp(model.matrix(mmmm)%*%mmmm$coef)+1))
 qresBoot[i,] = sort(sapply(1:(nrow(dados)), function(a) qnorm(pgeom(dddd$repl[a], mu[a]))))
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
     ylim = c(min(bands2[1,]),6),col=color,
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
plot(mu.hat,pearsonRes, pch = 16, cex = .6, xlab = expression(hat(mu)),
     ylab = "Residuo de Pearson", col = ifelse(pearsonRes > 3 |pearsonRes < -3, "red","black" ))
abline(h=3)
identify(mu.hat,pearsonRes)

# componentes do desvio
plot(mu.hat,devianceRes, pch = 16, cex = .6, xlab = expression(hat(mu)), 
     ylab = "Componentes do desvio", col = ifelse(devianceRes > 3 | devianceRes < -3, "red","black"),
     ylim = c(-3.1,3.1))
abline(h=c(-3,3))

# residuo quantilico
plot(mu.hat,qRes, pch = 16, cex = .6, xlab = expression(hat(mu)), 
     ylab = "Resíduo quantilico", col = ifelse(qRes > 3 | qRes < -3, "red", "black"),
     ylim = c(-3.1,3.1),type="n")
abline(h=c(-3,3))
points(mu.hat,qRes,pch = 16, cex = .6,col = ifelse(qRes > 3 | qRes < -3, "red", "black"))
identify(mu.hat,qRes)

# ***********************
# Correlacao serial (nao tenho certeza da ordem na real)
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
     ylab = "Resíduos de Pearson", xlab = expression(hat(eta)))

# Componentes do desvio
plot(design.X%*%coef.fit,devianceRes, pch = 16, cex = .6, 
     ylab = "Componentes do desvio", xlab = expression(hat(eta)))

# Quantilico
plot(design.X%*%coef.fit,qRes, pch = 16, cex = .6, 
     ylab = "Resíduo quantílico", xlab = expression(hat(eta)))

# ***************************
# Residuos x covariavel
# ***************************

# genero
plot(dados$gender,qRes, pch = 16, cex = .6, 
     ylab = "Resíduo quantílico", xlab = "Gênero", names = c("Feminino","Masculino"))

# nota de matematica
plot(dados$math,qRes, pch = 16, cex = .6, 
     ylab = "Resíduo quantílico", xlab = "Nota de matemática")

# programa
plot(dados$prog,qRes, pch = 16, cex = .6, 
     ylab = "Resíduo quantílico", xlab = "Programa")

# ---------------------------------------------------------------------------
# Pontos influentes
# ---------------------------------------------------------------------------

plot(mu.hat,diag(hat.matrix),pch=16,cex=.6,xlab=expression(hat(mu)),
     ylab="Alavancagem",ylim = c(.005,.05),
     col = ifelse(diag(hat.matrix)>2*p/n,"red","black"));abline(h=2*p/n)

plot(dffits(modelo),pch=16,cex=.6,col=ifelse(dffits(modelo)>2*sqrt(p/n),"red","black"),
     ylab = "DFFitts", xlab = "")
abline(h=2*sqrt(p/n))
identify(dffits(modelo))

cook.modelo = cooks.distance(modelo)
cook.f = pf(cook.modelo, df1 = p, df2 = n-p, lower.tail = F)
plot(cook.modelo, pch = 16, cex=.6, xlab ="",ylab="Distância de Cook",
     col = ifelse(cook.f <= .5, "red", "black"))
abline(h = qf(.5,df1=p,df2=n-p))

# ---------------------------------------------------------------------------
# Variacao parcial
# ---------------------------------------------------------------------------
alavanca = which(diag(hat.matrix)>2*p/n)
ajuste = which(dffits(modelo)>2*sqrt(p/n))

influentes = unique(c(alavanca,ajuste))


mod1 = glm(daysabs ~ math + gender + prog, data = dados[-influentes,],family = MASS::negative.binomial(theta = 1, link = "log"))
mod2 = glm(daysabs ~ math + gender + prog, data = dados[-influentes[-21],],family = MASS::negative.binomial(theta = 1, link = "log"))
mod3 = glm(daysabs ~ math + gender + prog, data = dados[-influentes[-22],],family = MASS::negative.binomial(theta = 1, link = "log"))
mod4 = glm(daysabs ~ math + gender + prog, data = dados[-influentes[-c(21,22)],],family = MASS::negative.binomial(theta = 1, link = "log"))

coef(mod1) # modelo sem
coef(mod2) # modelo com 1
coef(mod3) # modelo com 2
coef(mod4) # modelo com 1 e 2 
coef(modelo) # modelo com tudo

(coef(mod2)-coef(mod1))/coef(mod1)*100
(coef(mod3)-coef(mod1))/coef(mod1)*100
(coef(mod4)-coef(mod1))/coef(mod1)*100
(coef(modelo)-coef(mod1))/coef(mod1)*100
