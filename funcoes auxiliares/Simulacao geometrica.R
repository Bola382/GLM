rreggeom = function(covs, betas, invlink = exp){
 betas = as.vector(betas)
 covs = as.matrix(covs)
 n = nrow(covs)
 
 predlin = as.vector(covs%*%betas) # preditor linear
 mu = invlink(predlin) # media
 
 prob = 1/(1+mu) # probabilidade
 
 rgeom(n, prob)
}
