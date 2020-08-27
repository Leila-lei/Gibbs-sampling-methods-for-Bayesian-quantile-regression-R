
generate_data = function(n, d, dist_eps, beta_p){
  x=t(cbind(matrix(rep(1,n),n,1),matrix(rnorm(n*(d-1),0,1),n,d-1)))
  if (dist_eps=="stdN"){epsilon=matrix(rnorm(n,0,1),n,1)}
  else if (dist_eps=="heteroN"){epsilon=matrix(rnorm(n,0,1)*(1+x[3,]),n,1)}
  else if (dist_eps=="student"){epsilon=matrix(rt(n,3),n,1)}
  y = t(x)%*%beta_p + epsilon
  return(list(x, y))

}


