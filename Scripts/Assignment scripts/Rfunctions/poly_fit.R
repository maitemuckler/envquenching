# List the functions in this file
if(!exists('list_functions_polyfit')){
  print('Defining functions from poly_fit.R:')
  print('-------------------------------------')
  print(' 1. poly_fit(x, y, w, degree, nit, nsigma)')

  list_functions_polyfit = F
}

poly_fit <- function(x, y, w, degree, nit, nsigma){
  if(missing(degree)){degree = 1}
  if(missing(nit)){nit = 1}
  if(missing(nsigma)){nsigma = 1e5}
  if(missing(w)){w = y; w[] = 1}
  
  xx.xx = is.na(x) == F & is.na(y) == F
  x = x[xx.xx]
  y = y[xx.xx] 
  w = w[xx.xx]
  
  for(i in 1:nit){
    lm1 <- lm(y ~ poly(x, degree = degree, raw = T), weights = w)
    fit = lm1$fitted.values
    
    delta = y - fit
    chisqrt = sd(delta)
    
    x = x[abs(delta) < (nsigma * chisqrt)]
    y = y[abs(delta) < (nsigma * chisqrt)]
    w = w[abs(delta) < (nsigma * chisqrt)]
  }
  lm1 <- lm(y ~ poly(x, degree = degree, raw = T), weights = w)
  
  return(lm1)
}

