rregexp = function(covs, betas, invlink = exp){
 betas = as.vector(betas)
 covs = as.matrix(covs)
 n = nrow(covs)
 
 predlin = as.vector(covs%*%betas) # preditor linear
 mu = invlink(predlin) # media
 
 rate = 1/mu # taxa da exponencial
 
 rexp(n, rate = rate)
}