geoDev = function(resp,mu){
 n = length(resp)
 dev1 = function(mu) log(1+mu)
 dev2 = function(resp,mu){
  resp*log(resp) - resp*log(mu) - (resp+1)*log(1+resp) + (resp+1)*log(1+mu)
 }
 
 sum(sapply(1:n, function(a) ifelse(resp[a]==0, 2*dev1(mu[a]), 2*dev2(resp[a],mu[a]))))
}
