library(data.table)

### Definir diretórios ----
wdmain       <- "~/Work/Research/Astronomy/"
wdproject    <- paste0(wdmain, "Projects/environmental-quenching/")
wdcode       <- paste0(wdproject, "Scripts/")
wddata       <- paste0(wdmain, "Data/")

## Chamar outras funções ----
source(paste0(wdcode, 'Assignment/Rfunctions/poly_fit.R'))
source(paste0(wdcode, 'Assignment/Rfunctions/cosmodist.R'))

## Funções ----
angDist_rad <- function(ra1, dec1, ra2, dec2){
  if(length(ra1) >= length(ra2)){
    ra_v = ra1; dec_v = dec1
    ra0 = ra2; dec0 = dec2
  }else{
    ra_v = ra2; dec_v = dec2
    ra0 = ra1; dec0 = dec1    
  }
  xx.xx = ra_v == ra0 & dec_v == dec0
  temp = ra_v; temp[] = 0
  
  temp[!xx.xx] = acos(cos((90.0 - dec_v[!xx.xx]) * deg2rad) * cos((90.0 - dec0) * deg2rad) + 
                        sin((90.0 - dec_v[!xx.xx]) * deg2rad) * sin((90.0 - dec0) * deg2rad) * 
                        cos((ra_v[!xx.xx] - ra0) * deg2rad)) 
  
  return(temp)
}

## Constantes ----
deg2rad = 0.0174532925
rad2deg = 1 / 0.0174532925
c = 299792
arcsec2rad = deg2rad / 3600
H_0 = 70   # km s-1 Mpc-1
Omega_m = 0.3
Omega_l = 0.7

## Algumas opções ----
compute_dens = T
Rproj_max = 1 # Mpc
delta_v_max = 3000 #1500 #3000 used when comoving_distance = F --> maximum difference in velocity between the galaxy and its neighbour (in km/s)
lim_Mabs = -20.4

## Definir input e output data ----
input_file  <- "SDSS_DR18_Legacy_MGS_QSO"
output_file <- paste0('LocalDensity/DensityVolume/', input_file, '_density_Rproj', Rproj_max, '_Dv', delta_v_max, 'kms_Mabs', lim_Mabs, '.csv')

## Ler input data ----
input_data_i    <- fread(paste0(wddata, "Raw data/SDSS/", input_file, ".csv")) # targets (não ficam cortadas em mag: todas)
sdss_galaxies_i <- fread(paste0(wddata, "Raw data/SDSS/", input_file, ".csv")) # Galáxias que eu permito que sejam vizinhas (vou cortar em magAbs!)

## FIT DISTANCE ----
DM_Mpc_temp = vector()
zzs = seq(0.005, 0.2, 0.005)
i = 1
for(zz in zzs){
  temp <- cosmodist(zz, H0 = H_0, Omega_m = Omega_m, Omega_l = Omega_l)
  DM_Mpc_temp[i] = temp$DM_Mpc
  i = i + 1
}
DMz_fit <- poly_fit(zzs, DM_Mpc_temp, degree = 10, nit = 1, nsigma = 50)

sdss_galaxies_i$Dc = predict.lm(DMz_fit, newdata = data.frame(x = as.numeric(sdss_galaxies_i$z)))

## Cortes nos dados ----
xx.xx = sdss_galaxies_i$absPetro_r <= lim_Mabs 
sdss_galaxies = sdss_galaxies_i[xx.xx, ]

input_data = input_data_i#[1:10, ]

rm(input_data_i, sdss_galaxies_i)

print('##############################################################################')
print(sprintf('Galaxy catalogue redshift range: %.2f < z < %.2f', min(sdss_galaxies$z), max(sdss_galaxies$z)))
print('##############################################################################')

if(compute_dens){
  ###########################
  # COMPUTE DENSITY
  ###########################
  Nobj = nrow(input_data)
  igal = input_data$igal
  null_vector = igal
  null_vector[] = -99
  
  dens_galaxy = data.frame(igal, Ngals_vol = null_vector, logMstar_vol = null_vector, 
                           logdens_n = null_vector, logdens_mass = null_vector)
  
  progress = round(quantile(c(1:Nobj), probs = seq(0, 1, 0.05)))
  pp = 1
  
  id <- igal[1]
  
  for(id in igal){
    ra_i = input_data$ra[input_data$igal == id] # ra da galáxia com igal == id
    dec_i = input_data$dec[input_data$igal == id] # dec da galáxia com igal == id
    z_i = input_data$z[input_data$igal == id] # z da galáxia com igal == id

    Hz = H_0 * sqrt(Omega_m * (1 + z_i)**3 + Omega_l) # Hz da galáxia com z = z_i
    Dc_i = cosmodist(z_i, H0 = H_0, Omega_m = Omega_m, Omega_l = Omega_l)$DM_Mpc # Distância comóvel
    DA_i = cosmodist(z_i, H0 = H_0, Omega_m = Omega_m, Omega_l = Omega_l)$DA_Mpc # Angular diameter distance
    angDist_i = angDist_rad(ra_i, dec_i, sdss_galaxies$ra, sdss_galaxies$dec) # Distância angular entre a galáxia
                                                                              # e todas as outras da tabela
    
    dist_max_temp = Rproj_max / DA_i * rad2deg # degrees

    # Quais galáxias da tabela estão dentro desse critério?
    xx.xx_temp = angDist_i * rad2deg <= dist_max_temp & abs(sdss_galaxies$z - z_i) * c <= delta_v_max
    Ngals_vol = sum(xx.xx_temp) # Nr de galáxias que atendem o critério
    sdss_galaxies_temp = sdss_galaxies[xx.xx_temp, ] # Quais galáxias são essas? (inclui a própria galáxia analisada)
    logMstar_vol = log10(sum(10**sdss_galaxies_temp$lgm_tot_p50)) # Log da massa estelar dentro desse volume

    xx.xx = dens_galaxy$igal == id # posição da galáxia análisada na tabela dens_galaxy
    dens_galaxy$Ngals_vol[xx.xx] = Ngals_vol # adicionando quantas galáxias naquela volume
    dens_galaxy$logMstar_vol[xx.xx] = logMstar_vol # adicionando logmstar dentro do volume
    dens_galaxy$logdens_n[xx.xx] = log10((Ngals_vol + 1)/(pi*Rproj_max^2)) # Projetado?
    dens_galaxy$logdens_mass[xx.xx] = log10(10**logMstar_vol/(pi*Rproj_max^2))
    
    if(pp %in% progress){print(sprintf('%.0f / %.0f (%.1f%s)', pp, Nobj, pp/Nobj * 100, '%'))}
    pp = pp + 1
  }
  
  # temp_sdss = sdss_galaxies[, names(sdss_galaxies) %in% c('i', 'ra', 'dec', 'z')]

  #output = merge(input_data, dens_galaxy, by.x = 'igal', by.y = 'igal')
  output = dens_galaxy
  colnames(output)[-1] <- paste0(colnames(output)[-1], "_Rproj", Rproj_max)
  
  ###########################
  # WRITE OUTPUT
  ###########################
  write.csv(output, file = paste0(wddata, output_file), row.names = F, na = '')
}else{
  output = read.csv(paste0(wdinputdata, output_file), stringsAsFactors = F)
}
