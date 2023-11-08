# ========================================================================================
# Algoritmo escore de Fisher para a estimativa dos coefs de regressao
# com ligacao log
# ========================================================================================
library(emulator) # calcula formas quadraticas e outros produtos de vetor/matriz

# resp: vetor de variaveis resposta
# covs: matriz de covariaveis, para ter intercepto deve ter uns na 1ª coluna
# tol: valor do criterio de parada do algoritmo
# print: se == T, printa a iteracao e a diferenca ao valor anterior de beta

FSexp = function(resp,covs,tol=1e-8,print=T){
 n = length(resp)
 covs = as.matrix(covs)
 # ---------------
 # passo 0
 # ---------------
 # valor inicial de mu = resp
 mu = resp
 
 eta = log(mu) # preditor linear
 z = eta # resposta ajustada
 
 # ---------------------------
 # passo 1 do escore de Fisher
 # atualizando beta
 # ---------------------------
 beta = solve(crossprod(covs), crossprod(covs,z))
 
 iter = 1
 dif = Inf
 
 # -------------------------------------
 # passos 2 em diante do Score de Fisher
 # -------------------------------------
 while(dif > tol){
  eta = covs%*%beta[,iter]
  mu = exp(eta) 
  z = eta + (resp-mu)/mu
  
  beta = cbind(beta,solve(crossprod(covs), crossprod(covs,z)))
  iter = iter+1
  
  dif = sum(abs((beta[,iter]-beta[,iter-1])/(beta[,iter-1]+.1)))
  
  if(print==T) cat("Iteracao: ",iter, "\t Diferenca ao anterior: ", dif,"\r")
 }
 return(beta[,iter])
}