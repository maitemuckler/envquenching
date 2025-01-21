library(data.table)

### Definir diretórios ----
wdmain       <- "~/Work/Research/Astronomy/"
wdproject    <- paste0(wdmain, "Projects/environmental-quenching/")
wdcode       <- paste0(wdproject, "Scripts/")
wddata       <- paste0(wdmain, "Data/")
wdinputdata  <- paste0(wdproject, "Data/InputModel/")

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
find_Nth = T
delta_v_max = 3000 #1500  # used when comoving_distance = F --> maximum difference in velocity between the galaxy and its neighbour (in km/s)
Nth = 5 # Nth nearest neighbour
lim_Mabs = -20.4

## Definir input e output data ----
input_file  <- "SDSS_DR18_Legacy_MGS_QSO"
output_file <- paste0('LocalDensity/NthSpecGalaxy/', input_file, '_Nth_spec_galaxy_Neq', Nth, '_Dv', delta_v_max, 'kms_Mabs', lim_Mabs, '.csv')

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

print('##############################################################################')
print(sprintf('Galaxy catalogue redshift range: %.2f < z < %.2f', min(sdss_galaxies$z), max(sdss_galaxies$z)))
print('##############################################################################')

if(find_Nth){
  ###########################
  # FIND NTH NEAREST GALAXY
  ###########################
  Nobj = nrow(input_data)
  igal = input_data$igal
  null_vector = igal
  null_vector[] = -99
  
  Nth_galaxy = data.frame(igal, Nth_gal_igal = null_vector, Nth_gal_ra = null_vector, 
                          Nth_gal_dec = null_vector, Nth_gal_z = null_vector, 
                          Dist_comov = null_vector, Dist_proj = null_vector)
  
  progress = round(quantile(c(1:Nobj), probs = seq(0, 1, 0.05)))
  pp = 1
  
  for(id in igal){
    ra_i = input_data$ra[input_data$igal == id] # ra da galáxia com igal == id
    dec_i = input_data$dec[input_data$igal == id] # dec da galáxia com igal == id
    z_i = input_data$z[input_data$igal == id] # z da galáxia com igal == id

    Hz = H_0 * sqrt(Omega_m * (1 + z_i)**3 + Omega_l) # Hz da galáxia com z = z_i
    Dc_i = cosmodist(z_i, H0 = H_0, Omega_m = Omega_m, Omega_l = Omega_l)$DM_Mpc # Distância comóvel
    DA_i = cosmodist(z_i, H0 = H_0, Omega_m = Omega_m, Omega_l = Omega_l)$DA_Mpc # Angular diameter distance
    angDist_i = angDist_rad(ra_i, dec_i, sdss_galaxies$ra, sdss_galaxies$dec) # Distância angular entre a galáxia
                                                                              # e todas as outras da tabela
    dist_max_temp = 0.5 # degress
    repeat_find = T
    
    while(repeat_find == T & dist_max_temp <= 10){ # testar com 20
      angDist_i = angDist_rad(ra_i, dec_i, sdss_galaxies$ra, sdss_galaxies$dec)
      # xx.xx_temp = angDist_i* rad2deg * 3600 >= 5 & angDist_i * rad2deg <= dist_max_temp
      xx.xx_temp = sdss_galaxies$igal != id & angDist_i * rad2deg <= dist_max_temp # Critério 1: Galáxias que não seja a  
                                                                                   # galáxia analisada e que tenham 
                                                                                   # distância angular menor que
                                                                                   # dist_max_temp
      sdss_galaxies_temp = sdss_galaxies[xx.xx_temp, ] # subset só com as galáxias que atendem o critério
      angDist_i = angDist_i[xx.xx_temp] # distância angular das galáxias que atendem o critério
      dist_max_temp = 2 * dist_max_temp # atualiza a distância máxima
      
      mean_Dc = (sdss_galaxies_temp$Dc + Dc_i)/2 # Distância comóvel média entre a Dc de tds galáxias e da galáxia analisada
      Dist_comov = sqrt((sdss_galaxies_temp$Dc - Dc_i)**2 + (angDist_i * mean_Dc)**2) 
      xx.xx_vel = abs(sdss_galaxies_temp$z - z_i) * c <= delta_v_max # Critério 2: quais da tabela tem velocidade 
                                                                     # dentro do Critério 1?
      
      if(sum(xx.xx_vel) >= Nth){ # Se a quantidade de galáxias dentro dos Critérios 1 e 2 for maior ou igual ao número 
                                 # de vizinhos desejado
        repeat_find = F # Para de procurar
        Dist_proj = angDist_i * DA_i
        temp = sort(Dist_proj[xx.xx_vel], index.return = T)
        Nth_gal_id_i = sdss_galaxies_temp$igal[xx.xx_vel][temp$ix][Nth] # id da Nth galáxia
        
        #if((angDist_i[xx.xx_vel][temp$ix][1] * rad2deg * 3600) <= 5){Nth_gal_id_i = sdss_galaxies_temp$igal_lim[xx.xx_vel][temp$ix][Nth + 1]}
      }
    }
 
    if(sum(xx.xx_vel) >= Nth){
      Dist_proj = angDist_i * DA_i
      xx.xx_out = sdss_galaxies_temp$igal == Nth_gal_id_i # Posição da Nth galáxia na tabela sdss_galaxies_temp
      xx.xx = Nth_galaxy$igal == id # Posição da galáxia analisada na tabela Nth_galaxy
      Nth_galaxy$Nth_gal_igal[xx.xx] = sdss_galaxies_temp$igal[xx.xx_out] # Adicionando id da Nth galáxia em Nth_galaxy
      Nth_galaxy$Nth_gal_ra[xx.xx] = sdss_galaxies_temp$ra[xx.xx_out] # Adicionando ra da Nth galáxia em Nth_galaxy
      Nth_galaxy$Nth_gal_dec[xx.xx] = sdss_galaxies_temp$dec[xx.xx_out] # Adicionando dec da Nth galáxia em Nth_galaxy
      Nth_galaxy$Nth_gal_z[xx.xx] = sdss_galaxies_temp$z[xx.xx_out] # Adicionando z da Nth galáxia em Nth_galaxy
      Nth_galaxy$Dist_comov[xx.xx] = Dist_comov[xx.xx_out] # Adicionando Dist_comov da Nth galáxia em Nth_galaxy
      Nth_galaxy$Dist_proj[xx.xx] = (angDist_i * DA_i)[xx.xx_out] # Adicionando Dist_proj da Nth galáxia em Nth_galaxy
      Nth_galaxy$dens_comov[xx.xx] = ((Nth + 1)/((4/3)*pi*Dist_comov[xx.xx_out]^3)) # Densidade comóvel
      Nth_galaxy$dens_proj[xx.xx] = ((Nth + 1)/(pi*Dist_proj[xx.xx_out]^2)) # Densidade projetada
    }
    
    if(pp %in% progress){print(sprintf('%.0f / %.0f (%.1f%s)', pp, Nobj, pp/Nobj * 100, '%'))}
    pp = pp + 1
  }
  
  # temp_sdss = sdss_galaxies[, names(sdss_galaxies) %in% c('i', 'ra', 'dec', 'z')]

  #output = merge(input_data, Nth_galaxy, by.x = 'igal', by.y = 'igal')
  output = Nth_galaxy
  colnames(output)[-1] <- paste0(colnames(output)[-1], "_Neq", Nth)
  
  ###########################
  # WRITE OUTPUT
  ###########################
  write.csv(output, file = paste0(wddata, output_file), row.names = F, na = '')
}else{
  output = read.csv(paste0(wdinputdata, output_file), stringsAsFactors = F)
}
