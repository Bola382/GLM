pearsonResiduals = function(resp,mu){
 resp = as.vector(resp)
 mu = as.vector(mu)
 
 (resp-mu)/mu
}