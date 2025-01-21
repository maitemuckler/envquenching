Sigma_M_add <- function(Rhlr, e, logMstar){
  
  # área da elipse A = pi*a*b
  # df$Rhlr = a
  # df$e = 1 - (b/a) -> df$e = 1 - (b/df$Rhlr) -> ellipticity
  # Isolando para obter o b
  # df$e + (b/df$Rhlr) = 1
  # b = (1 - df$e) * df$Rhlr
  
  area    <- pi * (Rhlr) * ((1 - e) * Rhlr)
  Sigma_M <- (10^logMstar)/(area) #densidade superficial de massa (em unidades em M_Sun / kpc^2)
  return(Sigma_M)
  
}
