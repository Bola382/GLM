# ========================================================================================
# Algoritmo escore de Fisher para a estimativa dos coefs de regressao
# com ligacao log
# ========================================================================================
library(emulator) # calcula formas quadraticas e outros produtos de vetor/matriz

# resp: vetor de variaveis resposta
# covs: matriz de covariaveis, para ter intercepto deve ter uns na 1ª coluna
# eps: valor prox de zero a ser somado em y em caso de y=0 para podermos obter o valor
#      inicial de mu
# tol: valor do criterio de parada do algoritmo
# print: se == T, printa a iteracao e a diferenca ao valor anterior de beta

FSgeo = function(resp,covs,eps=1/6,tol=1e-8,print=T){
 n = length(resp)
 covs = as.matrix(covs)
 # ---------------
 # passo 0
 # ---------------
 mu = 1:n
 # valor inicial de mu = resp
 mu = sapply(mu, function(a) ifelse(resp[a]==0, resp[a] + eps, resp[a])) 
 
 eta = log(mu) # preditor linear
 peso = diag(mu/(1+mu)) # peso
 z = eta # resposta ajustada
 
 # ---------------------------
 # passo 1 do escore de Fisher
 # atualizando beta
 # ---------------------------
 beta = solve(quad.form(peso,covs), quad.3form(peso, covs, z))
 
 iter = 1
 dif = Inf
 
 # -------------------------------------
 # passos 2 em diante do Score de Fisher
 # -------------------------------------
 while(dif > tol){
  eta = covs%*%beta[,iter]
  mu = exp(eta) 
  diag(peso) = mu/(1+mu)
  z = eta + (resp-mu)/mu
  
  beta = cbind(beta,solve(quad.form(peso,covs), quad.3form(peso, covs, z)))
  iter = iter+1
  
  dif = sum(abs((beta[,iter]-beta[,iter-1])/(beta[,iter-1]+.1)))

  if(print==T) cat("Iteracao: ",iter, "\t Diferenca ao anterior: ", dif,"\r")
 }
 return(list(result = beta[,iter], historico = beta))
}