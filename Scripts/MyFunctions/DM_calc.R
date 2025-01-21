DM_calc <- function(redshifts){
  
  c       <- 299792
  H0      <- 70.2       # km/(s*Mpc)
  DH_Mpc  <- c / H0     # Mpc
  
  Omega_m <- 0.3
  Omega_l <- 0.7
  dz      <- 1e-6
  
  distance_modulus <- vector()
  
  for (i in 1:length(redshifts)) {
    
    z <- redshifts[i]
    
    print(i)
    
    zz  <- seq(0, z, dz)
    E_z <- sqrt(Omega_m * (1 + zz)**3 + Omega_l)
    
    DC_Mpc <- DH_Mpc * sum(dz / E_z) # (DC_Mpc = Comoving distance)
    DM_Mpc <- DC_Mpc       # Omega_z = 0 (flat Universe) (DM_Mpc = Transverse comoving distance)
    DL_Mpc <- (1 + z) * DM_Mpc # (DL_Mpc = luminosity distance)
    DM <- 5 * (log10(DL_Mpc * 1e6) - 1)
    
    distance_modulus <- append(distance_modulus, DM)
  }
  
  return(distance_modulus)
  
}