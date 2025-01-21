sigma_gal_correction <- function(sigma_gal, Rhlr, Scale){
  r_e <- Rhlr/Scale # fica em arcsec
  sigma_gal <- sigma_gal/((1.5/r_e)^(-0.066))
  return(sigma_gal)
}

# Simard+11
# SDSS structural parameters from n=4 bulge+disk decompositions (1123718 rows)