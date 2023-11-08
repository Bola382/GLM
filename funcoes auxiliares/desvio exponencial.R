expDev = function(resp,mu,comp=F){
 resp = as.vector(resp)
 mu = as.vector(mu)
 n = length(resp)
 if(comp==T){return(sign(resp-mu)*sqrt(2)*sqrt(resp/mu + log(mu) - log(resp)-1))}
 -2*n + 2*sum(resp/mu + log(mu) - log(resp))
}
