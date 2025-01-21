##################################################################################################
# PURPOSE:
#    Compute M_200 and R_200 (Delta relative to critical density) from 
#     M_180 and R_180 (Delta relative to mean density)
# CALLING SEQUENCE:
#    Delta180mto200c_NFW(m180_h, omega_m, H_0, z)
# INPUTS:
#    m180_h      input M_180 in units of h-1 Msun
#    omega_m     Omega_m (default = 0.238)
#    H_0         Hubble constant (default = 73)
#    z           redshift
# PARAMETERS:
#    
# OUTPUT:
#    logM180     log M_180 [M_sun]
#    r180        R_180 [Mpc]
#    logM200     log M_200 [M_sun]
#    r200        R_200 [Mpc]
#    c           concentration (Dutton & Maccio, 2014, MNRAS, 441, 3359)
# REQUIRED SCRIPTS:
#   ../../../../lib/Rlib/interp.R
##################################################################################################
#-------------------
# INTERP FUNCTION
#-------------------
#source('../../../../Rfunctions/interp.R')

#-------------------
# Delta180mto200c_NFW FUNCTION
#-------------------
Delta180mto200c_NFW <- function(m180_h, omega_m, H_0, z){

  if(missing(H_0)){H_0 = 73}                              # km s^-1 Mpc^-1
  if(missing(omega_m)){omega_m = 0.238}

  m200 = 10**seq(1, 20, 0.1)
  m180 = 10**m180_h / (H_0 / 100)                         # input m180 in units of h-1 Msun

  A = 10**(0.52 + (0.905 - 0.52) * exp(-0.617 * z**1.21)) # Dutton & Maccio, 2014, MNRAS, 441, 3359
  B = -0.101 + 0.026 * z                                  # Dutton & Maccio, 2014, MNRAS, 441, 3359

  G = 4.302e-3 * 1e-6                                     # Mpc * Msun^-1 * (km/s)^2

  m200_h = m200 * (H_0/100)
  c200 = A * (m200_h / 1e12)**B
  
  H_z = sqrt(H_0**2 * (omega_m * (1 + z)**3 + 1 - omega_m))
  omega_m_z = omega_m * (1 + z)**3 / (omega_m * (1 + z)**3 + 1 - omega_m)
  
  r180 = (m180 * 2 * G / (180 * omega_m * H_0**2))**(1/3) * (1 + z)**(-1)
  b200 = r180 * c200 * (m200 * 2 * G / (200 * H_z**2))**(-1/3)
  f200 = 180/2 * omega_m_z * (b200**3 / (log(1 + b200) - (b200 / (1 + b200))))
  
  fc200 = 200/2 * (c200**3 / (log(1 + c200) - (c200 / (1 + c200))))
  
  temp = (log10(fc200/f200))
  m200_out = interp(temp, 0, log10(m200)) 

  #plot(m200, f, log = 'xy', type = 'l', col = 'red')
  #lines(m200, fc)
  #abline(v = 10**m200_out)

  #plot(m200, temp, log = 'x')
  
  cout = A * (10**m200_out * (H_0/100) / 1e12)**B
  r200_out = (10**m200_out * 2 * G / (200 * H_z**2))**(1/3) 
  
  #data <- sprintf('%10s%10s%10s%10s%10s%10s', 'log M180', 'R180', 'log M200', 'R200', 'R180/R200', 'c')
  #print(data)

  #data <- sprintf('%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f', log10(m180), r180, m200_out, r200_out, r180/r200_out, cout)
  #print(data)

  # OUTPUT UNITS: MASSES --> Msun;  RADII --> Mpc
  output <- data.frame(logM180m = log10(m180), r180m = r180, 
                       logM200 = m200_out, r200 = r200_out, c200 = cout)
  return(output)
}

##################################################################################################
# PURPOSE:
#    Compute M_200 and R_200 (Delta relative to critical density) from 
#     M_200 and R_200 (Delta relative to mean density) --> NEW CATALOGUES YANG 
# CALLING SEQUENCE:
#    Delta180mto200c_NFW(m180_h, omega_m, H_0, z)
# INPUTS:
#    m200        input M_200_mean 
#    omega_m     Omega_m (default = 0.275)
#    H_0         Hubble constant (default = 70.2)
#    z           redshift
# PARAMETERS:
#    
# OUTPUT:
#    logM200m     log M_200_mean [M_sun]
#    r200m        R_200_mean [Mpc]
#    logM200c     log M_200_crit [M_sun]
#    r200c        R_200_crit [Mpc]
#    c200         concentration (Dutton & Maccio, 2014, MNRAS, 441, 3359)
# REQUIRED SCRIPTS:
#   ../../../../lib/Rlib/interp.R
##################################################################################################
#-------------------
# INTERP FUNCTION
#-------------------
#source('../../../../Rfunctions/interp.R')

#-------------------
# Delta180mto200c_NFW FUNCTION
#-------------------
Delta200mto200c_NFW <- function(m200m, omega_m, H_0, z){
  
  if(missing(H_0)){H_0 = 70.2}                            # km s^-1 Mpc^-1
  if(missing(omega_m)){omega_m = 0.275}
  
  m200c = 10**seq(1, 20, 0.1)
  m200m = 10**m200m                                       # input m200 relative to mean density
  
  A = 10**(0.52 + (0.905 - 0.52) * exp(-0.617 * z**1.21)) # Dutton & Maccio, 2014, MNRAS, 441, 3359
  B = -0.101 + 0.026 * z                                  # Dutton & Maccio, 2014, MNRAS, 441, 3359
  
  G = 4.302e-3 * 1e-6                                     # Mpc * Msun^-1 * (km/s)^2
  
  m200c_h = m200c * (H_0/100)
  c200 = A * (m200c_h / 1e12)**B
  
  H_z = sqrt(H_0**2 * (omega_m * (1 + z)**3 + 1 - omega_m))
  omega_m_z = omega_m * (1 + z)**3 / (omega_m * (1 + z)**3 + 1 - omega_m)
  
  r200m = (m200m * 2 * G / (180 * omega_m * H_0**2))**(1/3) * (1 + z)**(-1)
  b200 = r200m * c200 * (m200c * 2 * G / (200 * H_z**2))**(-1/3)
  f200 = 200/2 * omega_m_z * (b200**3 / (log(1 + b200) - (b200 / (1 + b200))))
  
  fc200 = 200/2 * (c200**3 / (log(1 + c200) - (c200 / (1 + c200))))
  
  temp = (log10(fc200/f200))
  m200c_out = interp(temp, 0, log10(m200c)) 
  
  #plot(m200, f, log = 'xy', type = 'l', col = 'red')
  #lines(m200, fc)
  #abline(v = 10**m200_out)
  
  #plot(m200, temp, log = 'x')
  
  cout = A * (10**m200c_out * (H_0/100) / 1e12)**B
  r200c_out = (10**m200c_out * 2 * G / (200 * H_z**2))**(1/3) 

  # OUTPUT UNITS: MASSES --> Msun;  RADII --> Mpc
  output <- data.frame(logM200m = log10(m200m), r200m = r200m, 
                       logM200c = m200c_out, r200c = r200c_out, c200 = cout)
  return(output)
}

##################################################################################################
# PURPOSE:
#    Compute M_100 and R_100 (Delta relative to critical density) from 
#     M_200 and R_200 (Delta also relative to critical density)
# CALLING SEQUENCE:
#    Delta200to100_NFW(logM200, omega_m, H_0, z)
# INPUTS:
#    logM200     input log M200 in units of log Msun
#    omega_m     Omega_m (default = 0.238)
#    H_0         Hubble constant (default = 73)
#    z           redshift (default = 0)
# PARAMETERS:
#    
# OUTPUT:
#    logM100     log M_100 [M_sun]
#    r100        r_100 [Mpc]
#    logM200     log M_200 [M_sun]
#    r200        r_200 [Mpc]
#    c100        concentration 
#    c200        concentration (Dutton & Maccio, 2014, MNRAS, 441, 3359)
# REQUIRED SCRIPTS:
#   ../../../../lib/Rlib/interp.R
##################################################################################################

Delta200to100_NFW <- function(logM200, omega_m, H_0, z){
  
  if(missing(H_0)){H_0 = 73}                              # km s^-1 Mpc^-1
  if(missing(omega_m)){omega_m = 0.238}
  if(missing(z)){z = 0}
  
  H_z = sqrt(H_0**2 * (omega_m * (1 + z)**3 + 1 - omega_m))
  omega_m_z = omega_m * (1 + z)**3 / (omega_m * (1 + z)**3 + 1 - omega_m)
  
  m100 = 10**seq(5, 16, 0.1)
  m200 = 10**logM200                                      # input in units of log Msun
  
  A = 10**(0.52 + (0.905 - 0.52) * exp(-0.617 * z**1.21)) # Dutton & Maccio, 2014, MNRAS, 441, 3359
  B = -0.101 + 0.026 * z                                  # Dutton & Maccio, 2014, MNRAS, 441, 3359
  
  G = 4.302e-3 * 1e-6                                     # Mpc * Msun^-1 * (km/s)^2
  
  m200_h = m200 * (H_0/100)
  c200 = A * (m200_h / 1e12)**B
  r200 = (m200 * 2 * G / (200 * H_z**2))**(1/3) 
  rs = r200 / c200

  r100 = (m100 * 2 * G / (100 * H_z**2))**(1/3) 
  c100 = r100 / rs
  f100 = 100 * (c100**3 / (log(1 + c100) - (c100 / (1 + c100))))
  
  f200 = 200 * (c200**3 / (log(1 + c200) - (c200 / (1 + c200))))

  temp = (log10(f100/f200))
  m100_out = interp(temp, 0, log10(m100)) 

  r100_out = (10**m100_out * 2 * G / (100 * H_z**2))**(1/3) 
  r200_out = (m200 * 2 * G / (200 * H_z**2))**(1/3) 
  
  c100_out = r100_out / rs
  c200_out = r200 / rs

  # OUTPUT UNITS: MASSES --> Msun;  RADII --> Mpc
  output <- data.frame(logM100 = m100_out, r100 = r100_out, 
                       logM200 = log10(m200), r200 = r200_out, c100 = c100_out, c200 = c200_out)
  return(output)
}

##################################################################################################
# PURPOSE:
#    Compute M_Delta2 and R_Delta2 (Delta relative to critical density) from 
#     M200 and R200 (Delta also relative to critical density)
# CALLING SEQUENCE:
#    Delta200toDelta2_NFW(logM200, Delta2, omega_m, H_0, z)
# INPUTS:
#    logM200     input log M200 in units of log Msun
#    Delta2      Output Delta relative to critical density
#    omega_m     Omega_m (default = 0.238)
#    H_0         Hubble constant (default = 73)
#    z           redshift (default = 0)
# PARAMETERS:
#    
# OUTPUT:
#    logM_Delta2 log M_Delta2 [M_sun]
#    r_Delta2    r_Delta2 [Mpc]
#    logM200     log M200 [M_sun]
#    r_200       r_200 [Mpc]
#    c_200       concentration (c200 from Dutton & Maccio, 2014, MNRAS, 441, 3359)
#    c_Delta2    concentration 
# REQUIRED SCRIPTS:
#   ../../../../lib/Rlib/interp.R
##################################################################################################

Delta200toDelta2_NFW <- function(logM200, Delta2, omega_m, H_0, z){
  
  if(missing(H_0)){H_0 = 73}                              # km s^-1 Mpc^-1
  if(missing(omega_m)){omega_m = 0.238}
  if(missing(z)){z = 0}
  
  H_z = sqrt(H_0**2 * (omega_m * (1 + z)**3 + 1 - omega_m))
  omega_m_z = omega_m * (1 + z)**3 / (omega_m * (1 + z)**3 + 1 - omega_m)
  
  mDelta2 = 10**seq(5, 16, 0.1)
  m200 = 10**logM200                                      # input in units of log Msun
  
  A = 10**(0.52 + (0.905 - 0.52) * exp(-0.617 * z**1.21)) # Dutton & Maccio, 2014, MNRAS, 441, 3359
  B = -0.101 + 0.026 * z                                  # Dutton & Maccio, 2014, MNRAS, 441, 3359
  
  G = 4.302e-3 * 1e-6                                     # Mpc * Msun^-1 * (km/s)^2
  
  m200_h = m200 * (H_0/100)
  c200 = A * (m200_h / 1e12)**B
  r200 = (m200 * 2 * G / (200 * H_z**2))**(1/3) 
  rs = r200 / c200
  
  rDelta2 = (mDelta2 * 2 * G / (Delta2 * H_z**2))**(1/3) 
  cDelta2 = rDelta2 / rs
  fDelta2 = Delta2 * (cDelta2**3 / (log(1 + cDelta2) - (cDelta2 / (1 + cDelta2))))
  
  f200 = 200 * (c200**3 / (log(1 + c200) - (c200 / (1 + c200))))
  
  temp = (log10(fDelta2/f200))
  mDelta2_out = interp(temp, 0, log10(mDelta2)) 
  
  rDelta2_out = (10**mDelta2_out * 2 * G / (Delta2 * H_z**2))**(1/3) 
  r200_out = (m200 * 2 * G / (200 * H_z**2))**(1/3) 
  
  cDelta2_out = rDelta2_out / rs
  c200_out = r200 / rs
  
  # OUTPUT UNITS: MASSES --> Msun;  RADII --> Mpc
  output <- data.frame(logmDelta2 = mDelta2_out, rDelta2 = rDelta2_out, 
                       logM200 = log10(m200), r200 = r200_out, cDelta2 = cDelta2_out, c200 = c200_out)
  return(output)
}
