# List the functions in this file
if(!exists('list_functions_cosmodist')){
  print('Defining functions in cosmodist.R:')
  print('-------------------------------------')
  print(' 1. cosmodist(z, H0, Omega_m, Omega_l, precision)')
  print(' 2. cosmodist_zList(z, H0, Omega_m, Omega_l, precision)')
  print('        --> Function 2 is not working yet!!!')
  
  list_functions_cosmodist = F
}

##################################################################################################
# FUNCTION cosmodist
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
cosmodist <- function(z, H0, Omega_m, Omega_l, dz){
  if(missing(dz)){dz = 1e-6} # incremento dz
  
  c      <- 299792     # velocidade da luz (km/s).
  DH_Mpc <- c / H0     # Hubble distance (Mpc).
  tH_Gyr <- 978 / H0   # Hubble time (Gyr).
  
  if(z > 0){
    
    zz  <- seq(0, z, dz)
    E_z <- sqrt(Omega_m * (1 + zz)**3 + Omega_l)
    
    DC_Mpc <- DH_Mpc * sum(dz / E_z) # Comoving distance (Mpc).
    DM_Mpc <- DC_Mpc                 # Proper motion distance para Omega_z = 0 (flat Universe) (Mpc).
    DA_Mpc <- DM_Mpc / (1 + z)       # Angular diameter distance (Mpc).
    DL_Mpc <- (1 + z) * DM_Mpc       # Luminosity distance (Mpc).
    
    TL_Gyr  <- tH_Gyr * sum(dz / ((1 + zz) * E_z)) # Lookback time (Gyr).
    DM      <- 5 * (log10(DL_Mpc * 1e6) - 1)       # Distance modulus
    VC_Gpc3 <- (4 * pi / 3) * DM_Mpc**3 * 1e-9     # Comoving volume para Ωk = 0 (Gpc³).
    
    cosmo = data.frame(DH_Mpc, tH_Gyr, DM_Mpc, DA_Mpc, DL_Mpc, VC_Gpc3, TL_Gyr, DM)
  }else{
    cosmo = data.frame(DH_Mpc = DH_Mpc, tH_Gyr = tH_Gyr, DM_Mpc = 0, DA_Mpc = 0, DL_Mpc = 0, VC_Gpc3 = 0, TL_Gyr = 0, DM = 0)
  }

  return(cosmo)
}


##################################################################################################
# FUNCTION cosmodist
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
cosmodist_zList_temp <- function(z, H0, Omega_m, Omega_l, dz_integral, z_min_fit, poly_order, dz_fit){
  source('poly_fit.R')
  if(missing(dz_integral)){dz_integral = 1e-6}
  if(missing(poly_order)){poly_order = 10}
  if(missing(z_min_fit)){z_min_fit = 0.005}
  if(missing(dz_fit)){dz_fit = 0.001}
  
  zz = seq(max(z_min_fit, min(z)), max(z), dz_fit)
  
  temp <- cosmodist(zz[1], H0, Omega_m, Omega_l, dz_integral)
  cosmo_temp = temp
  for(i in 2:length(zz)){
    temp <- cosmodist(zz[i], 73, 0.238, 0.762, dz_integral)
    cosmo_temp = rbind(cosmo_temp, temp)
  }
  
  DH_Mpc = c / H0     # Mpc
  tH_Gyr = 978 / H0   # Gyr
  
  yy = cosmo$DM_Mpc
  kk <- poly_fit(zz, yy, poly_order, 1, 50)
  tH_Gyr = predict.lm(kk, newdata = list(x = z))
  
  yy = cosmo$TL_Gyr
  kk <- poly_fit(zz, yy, poly_order, 1, 50)
  tH_Gyr = predict.lm(kk, newdata = list(x = z))
  
  DA_Mpc = DM_Mpc / (1 + z)
  DL_Mpc = (1 + z) * DM_Mpc
  
  DM = 5 * (log10(DL_Mpc * 1e6) - 1)
  
  VC_Gpc3 = 4 * pi / 3 * DM_Mpc**3 * 1e-9
  
  cosmo = data.frame(DH_Mpc, tH_Gyr, DM_Mpc, DA_Mpc, DL_Mpc, VC_Gpc3, TL_Gyr, DM)
  
  yy_i = c('DM_Mpc', 'DA_Mpc', 'DL_Mpc', 'VC_Gpc3', 'TL_Gyr', 'DM')
  
  cosmo[z <= z_min, names(cosmo) %in% yy_i] = 0
  
  return(cosmo)
}

##################################################################################################  
##################################################################################################  
cosmodist_zList <- function(z, H0, Omega_m, Omega_l, dz_integral, z_min_fit, poly_order, dz_fit){
  if(missing(dz_integral)){dz_integral = 1e-6}
  if(missing(poly_order)){poly_order = 10}
  if(missing(z_min_fit)){z_min_fit = 0.005}
  if(missing(dz_fit)){dz_fit = 0.001}
  
  dex = max(z) / max(z_min_fit, min(z))
  if(dex > 10){
    dex = log10(max(z) / max(z_min_fit, min(z)))
    
  }
  zz = seq(max(0.005, min(z)), max(z), 0.001)
  
  temp <- cosmodist(zz[1], H0, Omega_m, Omega_l, dz)
  cosmo_temp = temp
  for(i in 2:length(zz)){
    temp <- cosmodist(zz[i], 73, 0.238, 0.762)
    cosmo_temp = rbind(cosmo_temp, temp)
  }
  
  DH_Mpc = c / H0     # Mpc
  tH_Gyr = 978 / H0   # Gyr
  
  yy = cosmo$DM_Mpc
  kk <- poly_fit(zz, yy, poly_order, 1, 50)
  tH_Gyr = predict.lm(kk, newdata = list(x = z))
  
  yy = cosmo$TL_Gyr
  kk <- poly_fit(zz, yy, poly_order, 1, 50)
  tH_Gyr = predict.lm(kk, newdata = list(x = z))
  
  DA_Mpc = DM_Mpc / (1 + z)
  DL_Mpc = (1 + z) * DM_Mpc
  
  DM = 5 * (log10(DL_Mpc * 1e6) - 1)
  
  VC_Gpc3 = 4 * pi / 3 * DM_Mpc**3 * 1e-9
  
  cosmo = data.frame(DH_Mpc, tH_Gyr, DM_Mpc, DA_Mpc, DL_Mpc, VC_Gpc3, TL_Gyr, DM)
  
  yy_i = c('DM_Mpc', 'DA_Mpc', 'DL_Mpc', 'VC_Gpc3', 'TL_Gyr', 'DM')
  
  cosmo[z <= 0.005]
}


# Dear Marina,
# 
# Hera are my SM routines for luminosity distance.
# 
# cheers
# 
# Gary
# 
# P.S. Here is the relative precision of each approximation. I recommend Liu to be sure.
# 
# dlum 34	# cosmological luminosity distance (h-1 Mpc)
# # args: Omega_m^0 Omega_lambda^0 (scalars) z (vector) [Pen | WU | AK  approximation]
# # sources of approximations:
# # Pen 1999, ApJS 120, 49
# # Wickramasinghe & Ukwatta 2010, MNRAS 406, 548
# # Adachi & Kasai Prog. Theor. Phys. 127 (2012), 145-152, arXiv:1111.6396
# # Liu, Ma, Zhang & Yang 11, MNRAS 412, 2685
# # Liu is best with 10^-5 accuracy, but is 10x slower than others
# # then best is:
# # 0 < z < 0.14: Pen
# # 0.14 < z < 0.19: WU
# # 0.19 < z < 1.4: AK
# # 1.4 < z < 1.55: Pen
# # 1.55 < z < 5.8: AK
# # 5.8 < z < 7.6: Pen
# # z > 7.6; AK
# local set _zp1 = 1+$3
# if ($?4 && abs($1+$2-1) < 0.001) {
#   define coverH100 2997.9
#   if ('$4' == 'P' || '$4' == 'Pen' ) {  # Pen 99
#     set $0 = $coverH100*(1+$3)*(etaPen(1,$1)-etaPen(1/_zp1,$1))
#   } else { 
#     if ('$4' == 'WU') { # Wickramasinghe & Ukwatta 10
#       set $0 = $coverH100/3*_zp1/((1-$1)**(1/6)*$1**(1/3))*(PsiWU10(0,$1)-PsiWU10($3,$1))
#     } else {
#       if ('$4' == 'AK') { # Adachi & Kasai 12
#         set $0 = 2*$coverH100*_zp1/sqrt($1)*(phiAK($1,0)-phiAK($1,$3)/sqrt(_zp1))
#       } else {
#         if ('$4' == 'Liu') { # Liu et al. 11
#           local define _s (((1-$1)/$1)**(1/3))
#           set $0 = $coverH100*_zp1/sqrt($_s*$1)*(T_Liu($_s)-T_Liu($_s/_zp1))
#         } else {
#           echo DLUM: cannot recognize arg-4 = $4
#           return
#         }
#       }
#     }
#   }
# } else {
#   set $0 = dpm($1,$2,$3)*(1+$3)
# }
# 
# phiAK 2 # kernel for luminosity distance approximation by Adachi & Kasai (2012)
# # args: Omega_m z
# foreach vec (x x2 x3) {
#   set _$vec local
# }
# set _x = (1-$1)/$1/(1+$2)**3
# set _x2 = _x*_x
# set _x3 = _x2*_x
# set $0 = (1+1.320*_x+0.4415*_x2+0.02656*_x3)/(1+1.392*_x+0.5121*_x2+0.03944*_x3)
# 
# T_Liu 1 # kernel for luminosity distance calculation by Liu et al. (2011)
# local set _m = 2*sqrt($1*$1-$1+1)/$1 + 2/$1 - 1
# set $0 = 4*R_F(_m, _m+3+sqrt(12), _m+3-sqrt(12))
# 
# overdensityfromb 2 # overdensity (relative to Universe) given FoF b
# # args: b_FoF c=r_Delta/r_s
# # source: More et al. 2011, ApJS 195, 4, eqs. (13) & (14)
# set $0 = 244.86/(5*$1)**3*(ln($2+1)-$2/($2+1))*(1+$2)**2/$2**2 - 1
# 
# etaPen 2 # auxiliary function for luminosity distances by Ue-Li Pen
# # args: a=1/(1+z), Omega_m^0
# local set s3 = (1-$2)/$2
# local set s = s3**(1/3)
# set $0 = 2*sqrt(s3+1)*(1/$1**4-0.154*s/$1**3+0.4304*s*s/$1**2+0.19097*s3/$1+0.066941*s3*s)**(-1/8) 
# 
# PsiWU10 2 # auxiliairy function for luminosity distances by Wickramasinghe & Ukwatta 10
# # args: z Omega_m^0
# local set alpha = 1 + 2*(1-$2)/$2/(1+$1)**3
# local set x = ln(alpha+sqrt(alpha*alpha-1))
# set $0 = 3*2**(2/3)*x**(1/3)*(1-x*x/252+x**4/21060)

