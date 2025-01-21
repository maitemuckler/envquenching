
# List the functions in this file
if(!exists('list_functions_NFW')){
  print('Defining functions from NFW.R:')
  print('-------------------------------------')
  print(' 1. getN_NFW(c_vir, R)')
  print(' 2. getSigma_NFW(R, r_vir, c_vir, Mv)')
  print(' 3. getRho3D_NFW(r, r_vir, c_vir, Mv)')
  print(' 4. getSigma3D_NFW(r, r_vir, c_vir, Mv)')
  print(' 5. getSigmaSphere_NFW(R, rmax, r_vir, c_vir, Mv)')
  print(' 6. getNSphere_NFW(R, c_vir, r_vir, rmax)')
  print(' 7. getCorr_factor(c_vir, breaks, mids, dbin)')
  print(' 8. Sigma_NFW(c_vir, R, Rmin, Rmax, p_maggie, ddd, rho, SigmaSphere, r_max)')
  print(' 9. getRho3D_Ein(r, r_vir, c_vir, Mv)')
  print('10. Sigma_Ein(c_vir, n, R, Rmin, Rmax, p_maggie)')
  
  list_functions_NFW = F
}

##################################################################################################
# 1. FUNCTION getN_NFW
##################################################################################################
# PURPOSE:
#    
# CALLING SEQUENCE:
#    
# INPUTS:
#    
# PARAMETERS:
#    
# OUTPUT:
#    
# REQUIRED SCRIPTS:
#   
##################################################################################################  
getN_NFW <- function(c_vir, R){
  r_vir = 1
  rs = r_vir / c_vir
  Rtilda = R / r_vir
  gc = 1 / (log(1 + c_vir) - c_vir / (1 + c_vir))
  x = 1 / (c_vir * Rtilda)
  Cm1 = x ; Cm1[] = 0
  Cm1[R > rs] = acos(x[R > rs])
  Cm1[R < rs] = acosh(x[R < rs])
  
  cc = gc * (Cm1 / sqrt(abs(c_vir**2 * Rtilda**2 - 1)) + log(c_vir * Rtilda/2))
  cc[c_vir == 1] = gc[c_vir == 1] * (1 + log(1 / 2))
  N = as.numeric(cc)
  return(N)
}

##################################################################################################
# 2. FUNCTION getSigma_NFW
##################################################################################################
# PURPOSE:
#    
# CALLING SEQUENCE:
#    
# INPUTS:
#    
# PARAMETERS:
#    
# OUTPUT:
#    
# REQUIRED SCRIPTS:
#   
##################################################################################################  
getSigma_NFW <- function(R, r_vir, c_vir, Mv){
  rs = r_vir / c_vir
  Rtilda = R / r_vir
  gc = 1 / (log(1 + c_vir) - c_vir / (1 + c_vir))
  x = 1 / (c_vir * Rtilda)
  Cm1 = x ; Cm1[] = 0
  Cm1[R > rs] = acos(x[R > rs])
  Cm1[R < rs] = acosh(x[R < rs])
  
  Sigma = (c_vir**2 * gc / (2 * pi)) * (Mv / r_vir**2) *
    ((1 - (abs(c_vir**2 * Rtilda**2 - 1))**(-1/2) * Cm1) / (c_vir**2 * Rtilda**2 - 1)) 
  
  Sigma[R == rs] = (c_vir**2 * gc * Mv / (6 * pi * r_vir**2))
  return(Sigma)
}

##################################################################################################
# 3. FUNCTION getRho3D_NFW
##################################################################################################
# PURPOSE:
#    
# CALLING SEQUENCE:
#    
# INPUTS:
#    
# PARAMETERS:
#    
# OUTPUT:
#    
# REQUIRED SCRIPTS:
#   
##################################################################################################  
getRho3D_NFW <- function(r, r_vir, c_vir, Mv){
  rs = r_vir / c_vir
  rtilda = r / r_vir
  gc = 1 / (log(1 + c_vir) - c_vir / (1 + c_vir))
  
  rho = Mv / r_vir**3 * c_vir**3 / (4 * pi) * gc / ((r / rs) * (1 + r / rs)**2)

  return(rho)
}

##################################################################################################
# 4. FUNCTION getSigma3D_NFW
##################################################################################################
# PURPOSE:
#    
# CALLING SEQUENCE:
#    
# INPUTS:
#    
# PARAMETERS:
#    
# OUTPUT:
#    
# REQUIRED SCRIPTS:
#   
##################################################################################################  
getSigma3D_NFW <- function(r, r_vir, c_vir, Mv){
  rs = r_vir / c_vir
  rtilda = r / r_vir
  gc = 1 / (log(1 + c_vir) - c_vir / (1 + c_vir))
  
  rho = Mv / r_vir**3 * c_vir**3 / (4 * pi) * gc / ((r / rs) * (1 + r / rs)**2)
  # Sigma  = (4 pi r^2 dr) * rho(r) / (2 pi r dr) = 2 r rho(r)
  Sigma = 2 * r * rho
  
  return(Sigma)
}

##################################################################################################
# 5. FUNCTION getSigmaSphere_NFW (Mamon+ 2010, A&A, 520, 30)
##################################################################################################
# PURPOSE:
#    
# CALLING SEQUENCE:
#    
# INPUTS:
#    
# PARAMETERS:
#    
# OUTPUT:
#    
# REQUIRED SCRIPTS:
#   
##################################################################################################  
getSigmaSphere_NFW <- function(R, rmax, r_vir, c_vir, Mv){
  if(missing(rmax)){rmax = r_vir}
  rs = r_vir / c_vir
  cmax = rmax / rs
  gc = 1 / (log(1 + c_vir) - c_vir / (1 + c_vir))
  X = R * c_vir / r_vir
  
  xx.xx1 = X > 0 & X < 1
  xx.xx2 = X == 1 & X < cmax
  xx.xx3 = X > 1 & X < cmax
  xx.xx4 = X == 0 | X > cmax
  X1 = X[xx.xx1]; X2 = X[xx.xx2]; X3 = X[xx.xx3]
  
  Stilda = X ; Stilda[] = 0
  Stilda[xx.xx1] = ((1 / ((1 - X1**2)**(3/2))) * acosh((cmax + X1**2) / ((cmax + 1) * X1)) - (1 / (cmax + 1)) * 
                      (sqrt(cmax**2 - X1**2) / (1 - X1**2)))
  
  Stilda[xx.xx2] = (((sqrt(cmax**2 - 1) * (cmax + 2)) / (3 * (cmax + 1)**2)) + 
                      (((-2 * cmax**3 - 4 * cmax**2 - cmax + 2) * (X2 - 1)) / (5 * (cmax + 1)**2 * sqrt(cmax**2 - 1))))
  
  Stilda[xx.xx3] = ((1 / (cmax + 1)) * ((sqrt(cmax**2 - X3**2)) / (X3**2 - 1)) - (1 / ((X3**2 - 1)**(3/2))) * acos((cmax + X3**2) / ((cmax + 1) * X3)))
  
  Stilda[xx.xx4] = 0
  
  Stilda = (1 / (2 * log(2) - 1)) * Stilda
  
  Sigma = (c_vir**2 / (pi * r_vir**2)) * gc * (log(2) - 1/2) * Mv * Stilda
  
  return(Sigma)
}

##################################################################################################
# 6. FUNCTION getNSphere_NFW (Mamon+ 2010, A&A, 520, 30)
##################################################################################################
# PURPOSE:
#    
# CALLING SEQUENCE:
#    
# INPUTS:
#    
# PARAMETERS:
#    
# OUTPUT:
#    
# REQUIRED SCRIPTS:
#   
##################################################################################################  
getNSphere_NFW <- function(R, c_vir, r_vir, rmax){
  if(missing(rmax)){rmax = r_vir}
  rs = r_vir / c_vir
  cmax = rmax / rs
  gc = 1 / (log(1 + c_vir) - c_vir / (1 + c_vir))
  X = R * c_vir / r_vir

  xx.xx1 = X == 0
  xx.xx2 = X > 0 & X < 1 & X < cmax
  xx.xx3 = X > 1 & X < cmax
  xx.xx4 = X == 1 & X < cmax
  xx.xx5 = X >= cmax
  X1 = X[xx.xx1]; X2 = X[xx.xx2]; X3 = X[xx.xx3]; X4 = X[xx.xx4]; X5 = X[xx.xx5]
  Mtilda = X ; Mtilda[] = 0

  if(sum(xx.xx1) > 0){Mtilda[xx.xx1] = 0}
  
  if(sum(xx.xx2) > 0){Mtilda[xx.xx2] = (sqrt(cmax**2 - X2**2) - cmax) / (cmax + 1) + log(((cmax + 1) * (cmax - sqrt(cmax**2 - X2**2))) / X2) + 
    (1 / sqrt(1 - X2**2)) * acosh((cmax + X2**2) / (X2 * (cmax + 1)))}
  
  if(sum(xx.xx3) > 0){Mtilda[xx.xx3] = (sqrt(cmax**2 - X3**2) - cmax) / (cmax + 1) + log(((cmax + 1) * (cmax - sqrt(cmax**2 - X3**2))) / X3) + 
    (1 / sqrt(X3**2 - 1)) * acos((cmax + X3**2) / (X3 * (cmax + 1)))}
  
  if(sum(xx.xx4) > 0){Mtilda[xx.xx4] = log((cmax + 1) * (c - sqrt(cmax**2 - 1))) - cmax / (cmax + 1) + 2 * sqrt((cmax - 1) / (cmax + 1))}
  
  if(sum(xx.xx5) > 0){Mtilda[xx.xx5] = log(cmax + 1) - cmax / (cmax + 1)}
  
  Mtilda = (1 / (log(2) - 1/2)) * Mtilda
  
  Msphere = gc * (log(2) - 1/2) * Mtilda
  
  return(Msphere)
}

##################################################################################################
# 7. FUNCTION getCorr_factor
##################################################################################################
# PURPOSE:
#    
# CALLING SEQUENCE:
#    
# INPUTS:
#    
# PARAMETERS:
#    
# OUTPUT:
#    
# REQUIRED SCRIPTS:
#   
##################################################################################################  
getCorr_factor <- function(c_vir, breaks, mids, dbin){
  int_f_r = vector()
  for(i in 1:length(mids)){
    #ddR = 0.0001
    #RR = seq(breaks[i], breaks[i + 1], ddR)
    #int_f_r[i] = sum(RR * getSigma_NFW(RR, 1, c_vir, 1)) * ddR
    int_f_r[i] = integrate(function(RR){RR * getSigma_NFW(RR, 1, c_vir, 1)}, lower = breaks[i], upper = breaks[i + 1])$value
  }
  f_r = getSigma_NFW(mids, 1, c_vir, 1) * (mids)**2 * dbin * log(10)
  gamma = (int_f_r / f_r)
  return(gamma)
}

##################################################################################################
# 8. FUNCTION Sigma_NFW
##################################################################################################
# PURPOSE:
#    
# CALLING SEQUENCE:
#    
# INPUTS:
#    
# PARAMETERS:
#    
# OUTPUT:
#    
# REQUIRED SCRIPTS:
#   
##################################################################################################  
Sigma_NFW <- function(c_vir, R, Rmin, Rmax, p_maggie, ddd, rho, SigmaSphere, r_max){
  if(missing(SigmaSphere)){SigmaSphere = F}
  if(missing(r_max)){r_max = 1}
  if(missing(p_maggie)){p_maggie = 1; SigmaSphere = F}else{SigmaSphere = T}
  if(c_vir > 0.1 & c_vir < 10){
    if(missing(ddd)){ddd = F}
    if(missing(rho)){rho = F}
    R = R[R >= Rmin & R <= Rmax]
    
    #dRR = 0.0001
    #RR = seq(Rmin, Rmax, dRR)
    if(SigmaSphere == F){
      if(ddd == F){
        #int_P = sum(2 * pi * RR * getSigma_NFW(RR, 1, c_vir, 1)) * dRR
        int_P = integrate(function(RR){2 * pi * RR * getSigma_NFW(RR, 1, c_vir, 1)}, lower = Rmin, upper = Rmax)$value
        pp <- getSigma_NFW(R, 1, c_vir, 1) 
      }else{
        if(rho == F){
          #int_P = sum(2 * pi * RR * getSigma3D_NFW(RR, 1, c_vir, 1)) * dRR
          int_P = integrate(function(RR){2 * pi * RR * getSigma3D_NFW(RR, 1, c_vir, 1)}, lower = Rmin, upper = Rmax)$value
          pp <- getSigma3D_NFW(R, 1, c_vir, 1)
        }else{
          # 3D DENSITY RHO
          #int_P = sum(4 * pi * RR**2 * getRho3D_NFW(RR, 1, c_vir, 1)) * dRR
          int_P = integrate(function(RR){4 * pi * RR**2 * getRho3D_NFW(RR, 1, c_vir, 1)}, lower = Rmin, upper = Rmax)$value
          pp <- getRho3D_NFW(R, 1, c_vir, 1) 
        }
      }
    }else{
      # NFW TRUNCATED AT R_MAX R_VIR 
      #int_P = sum(2 * pi * RR * getSigmaSphere_NFW(RR, r_max, 1, c_vir, 1)) * dRR
      int_P = integrate(function(RR){2 * pi * RR * getSigmaSphere_NFW(RR, r_max, 1, c_vir, 1)}, lower = Rmin, upper = Rmax)$value
      pp <- getSigmaSphere_NFW(R, 1, c_vir, 1) 
    }
    if(rho == F){
      p_r = 2 * pi * R * pp / int_P 
    }else{
      p_r = 4 * pi * R**2 * pp / int_P 
    }
    
    lnlike = -sum((p_maggie * log(p_r))[!is.na(log(p_r))])
  }else{
    lnlike = 1e50
  }
  #if(is.na(lnlike)){lnlike = 1e50}
  return(lnlike)
}


##################################################################################################
# 9. FUNCTION getRho3D_Ein
##################################################################################################
# PURPOSE:
#    
# CALLING SEQUENCE:
#    
# INPUTS:
#    
# PARAMETERS:
#    
# OUTPUT:
#    
# REQUIRED SCRIPTS:
#   
##################################################################################################  
getRho3D_Ein <- function(r, r_vir, c_vir, Mv, n){
  rtilda = r / r_vir
  
  # Here are the Einasto formulae in virial units.
  #
  # rho(r) = M_v/r_v^3 rho^hat(r/r_v,c,n)
  # rho^hat(x,c,n) = (2n)^(3n–1) c^3 / [2 pi gamma(3n,2n c^{1/n})] Exp[–2n (c x)^{1/n}] 
  #
  # where gamma is the incomplete gamma function (from 0 to x).
  #  
  # M(r) = M_vir M^hat(r/r_v,c,n)
  # M^hat(x,c,n) = gamma[3n,2n (cx)^{1/n}] / gamma[3n,2n c^{1/n}]
  
  aa = 3 * n; xx = 2 * n * c_vir^(1/n)
  gg = pgamma(xx, aa) * gamma(aa)
  rho_hat = (2 * n)**(3 * n - 1) * c_vir**3 / (2 * pi * gg) * exp(-2 * n * (c_vir * rtilda)**(1/n))
  rho = Mv / r_vir**3 * rho_hat 

  return(rho)
}

##################################################################################################
# 10. FUNCTION Sigma_NFW
##################################################################################################
# PURPOSE:
#    
# CALLING SEQUENCE:
#    
# INPUTS:
#    
# PARAMETERS:
#    
# OUTPUT:
#    
# REQUIRED SCRIPTS:
#   
##################################################################################################  
Sigma_Ein <- function(c_vir, n, R, Rmin, Rmax, p_maggie){
  #if(c_vir > 0.1 & c_vir < 10 & n > 0.1 & n < 10){
    if(missing(p_maggie)){p_maggie = 1; SigmaSphere = F}else{SigmaSphere = T}
    R = R[R >= Rmin & R <= Rmax]
    
    #dRR = 0.0001
    #RR = seq(Rmin, Rmax, dRR)
    
    # 3D DENSITY RHO
    #int_P = sum(4 * pi * RR**2 * getRho3D_Ein(RR, 1, c_vir, 1, n)) * dRR
    int_P = integrate(function(RR){4 * pi * RR**2 * getRho3D_Ein(RR, 1, c_vir, 1, n)}, lower = Rmin, upper = Rmax)$value
    pp <- getRho3D_Ein(R, 1, c_vir, 1, n) 

    p_r = 4 * pi * R**2 * pp / int_P 
    
    lnlike = -sum(p_maggie * log(p_r))
    #if(is.na(lnlike)){lnlike = 1e50}
  #}else{
  #  lnlike = 1e50
  #}
  return(lnlike)
}
