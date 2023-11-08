rregnbinom = function(covs, betas, phi, invlink = exp){
 betas = as.vector(betas)
 covs = as.matrix(covs)
 n = nrow(covs)
 
 predlin = as.vector(covs%*%betas) # preditor linear
 mu = invlink(predlin) # media
 
 rnbinom(n, size = phi, mu = mu)
}