# ===================================================================================
# Ajustando um MLG com resposta geometrica
# ===================================================================================
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())

data("ships", package = "MASS")
dados = ships

head(dados)

str(dados)

dados$year = factor(dados$year)
dados$period = factor(dados$period)

str(dados)

# observacao importante: sao varios navios daquele tipo, construidos em certo ano,
#                        com operacao em certo periodo por linha
#                        entao a variavel service conta os meses totais de operacao dos navios
#                        por isso os valores extremamente altos
#type: tipo de navio (algum dos navios e mais fragil?)
#year: ano de construcao (construcao piorou em algum ano?)
#      1960-64: 60
#      1965-69: 65
#      1970-74: 70
#      1975-79: 75
#period: periodo de operacao (imagino que seja a tripulacao dos navios)
#service: quantos meses ficaram em servico (obviamente mais tempo implica em mais acidentes)
#incidents: numero de acidentes com dano

summary(dados)
# importante lembrar que estamos considerando independencia entre os navios
# entao estamos assumindo que a operacao em um periodo nao afeta a operacao num outro

# Fonte: 
# P. McCullagh and J. A. Nelder, (1983), Generalized Linear Models (2nd). Chapman & Hall, section 6.3.2, page 205

# -----------------------------------------------------------------------
# Preprocessamento
# -----------------------------------------------------------------------

dados$service = dados$service/12 # agora sao anos agregados de servico

dados[which(dados$service==0)[-5],] # incidents = service = 0 pois o navio ainda nao tinha sido construido
dados[34,] # incidents = services = 0 por acidente, na verdade deve ser um NA (dito na fonte dos dados)

dados[34,"service"] = NA; dados[34,"incidents"] = NA

summary(dados)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Imputando dados faltantes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(mice) # multiple inputation with chained equations

init = mice(dados, maxit=0) 
meth = init$method # predictive mean matching
predM = init$predictorMatrix

set.seed(1)
imputed = mice(dados, method = meth, predictorMatrix = predM, m = 10)
dados.input = complete(imputed)

summary(dados.input)

# -------------------------------------------------------------------------
# Modelo
# -------------------------------------------------------------------------

modelo = glm(incidents ~ service + (type+year+period)^2, data = dados.input,family = MASS::negative.binomial(theta = 1, link = "log"))
design.X = model.matrix(modelo)
coef.fit = coef(modelo)

source("Funcoes auxiliares/escore de Fisher MLG geometrica.R")
modelo.meu = FSgeo(dados.input$incidents, design.X)

cbind(coef.fit,modelo.meu$result)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Redução de covariaveis
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# correlacao entre variaveis
boxplot(dados.input$service~dados.input$period)
boxplot(dados.input$service~dados.input$year)
boxplot(dados.input$service~dados.input$type)
table(dados.input$period,dados.input$year)
table(dados.input$period,dados.input$type)
table(dados.input$period,dados.input$year)
table(dados.input$year,dados.input$type)

# selecao de variaveis
stepAIC(modelo, direction = "both")

modelo = glm(incidents ~ service + type + year + period + type:year + 
              year:period, data = dados.input,family = MASS::negative.binomial(theta = 1, link = "log"))
design.X = model.matrix(modelo)
coef.fit = coef(modelo)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Diagnostico
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n = nrow(dados.input) # amostras
p = length(coef.fit) # numero de parametros

# estimativa da media g^(-1)(eta.hat)
mu.hat = exp(design.X%*%coef.fit)

# Estatistica de Pearson
sum((dados.input$incidents-mu.hat)^2/(mu.hat+mu.hat^2))
sum(residuals(modelo, type = "pearson")^2)

qchisq(0.05, df = n-p, lower.tail = F) # ok

sum(residuals(modelo, type = "pearson")^2)/(n-p) # estimativa de a(phi)

# Deviance
source("Funcoes auxiliares/desvio geometrica.R")
geoDev(dados.input$incidents,mu.hat)
deviance(modelo)
sum(residuals(modelo, type = "deviance")^2)

qchisq(0.05, df = n-p, lower.tail = F) # ok

deviance(modelo)/(n-p) # estimativa de a(phi)

# **************************
# Residuos
# **************************

# de Pearson
pearsonRes = residuals(modelo, type="pearson")

# Deviance
devianceRes = residuals(modelo, type = "deviance")

# Quantilicos
p.hat = c(1/(mu.hat+1))
qRes = sapply(1:(nrow(dados)), function(a) qnorm(pgeom(dados.input$incidents[a], p.hat[a])))

# *******************
# envelope
# *******************
source("Funcoes auxiliares/Simulacao geometrica.R")

B = 1000
qresBoot = matrix(NA, nrow = B, ncol = nrow(dados))
pearsonResBoot = devianceResBoot = qresBoot 

set.seed(2)
for(i in 1:B){
 repl = rreggeom(design.X,coef.fit)
 dddd = cbind.data.frame(repl,design.X[,-1])
 mmmm = glm(repl ~ ., data = dddd,maxit = 300,family = MASS::negative.binomial(theta = 1, link = "log"))
 
 mu = c(1/(exp(model.matrix(mmmm)%*%mmmm$coef)+1))
 qresBoot[i,] = sort(sapply(1:(nrow(dados)), function(a) qnorm(pgeom(dddd$repl[a], mu[a]))))
 pearsonResBoot[i,] = sort(residuals(mmmm,type = "pearson"))
 devianceResBoot[i,] = sort(residuals(mmmm,type = "deviance"))
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
     #ylim = c(min(bands2[1,]),max(bands[2,])),
     col=color,
     pch = 16, cex=size)
lines(qnorm((1:n-.5)/n) ,bands[1,], col = "gray50")
lines(qnorm((1:n-.5)/n) ,bands[2,], col = "gray50")
lines(qnorm((1:n-.5)/n), colMeans(qresBoot), lty = 2)
lines(qnorm((1:n-.5)/n) ,bands2[1,])
lines(qnorm((1:(nrow(dados))-.5)/nrow(dados)) ,bands2[2,])
legend("topleft", legend = c("Range", "IC 95%", "Média"), title = "Resíduo simulado", lty = c(1,1,2), col = c("black", "gray50", "black"))

# dist parece ok a nao ser por 1 ponto no qqplot dos residuos quantilicos
which.min(qRes) #observacao 22
dados.input[22,]


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

