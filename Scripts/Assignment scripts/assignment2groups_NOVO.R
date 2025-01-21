### Bibliotecas ----
library(data.table)
library(dplyr)
library(polynom)

### Definir variáveis ----
Rlim            <-  2.5   
absMag_lim      <- -12  
N_rvir          <-  20    
Ngals_min_group <-  1      

zmax <- 0.07
zmin <- 0.01
Mhalo_min_group <- 14

### Definir quantidades cosmológicas ----
DELTA <- 100
eta   <- 0.65   # eta = sigma / v_vir

deg2rad    <- 0.0174532925
rad2deg    <- 1 / 0.0174532925
H_0        <- 70.2   # km s-1 Mpc-1
G_cosmo    <- 4.302e-9 # Mpc Msun-1 (km s-1)^2 # DOIS G?
Omega_m    <- 0.3
Omega_l    <- 0.7
rho_mean_0 <- Omega_m * 3 * H_0**2 / (8 * pi * G_cosmo) # Msun Mpc-3

G       <- 6.67408e-11 # m3 kg-1 s-2 # DOIS G?
Msun2kg <- 1.98855e30 # kg
Mpc2m   <- 3.086e22   # m

### Definir diretórios ----
wdmain       <- "~/Work/Research/Astronomy/"
wdproject    <- paste0(wdmain, "Projects/environmental-quenching/")
wdcode       <- paste0(wdproject, "Scripts/")
wddata       <- paste0(wdmain, "Data/")
wdassigndata <- paste0(wddata, "environmental-quenching-data/Assignment2groups/")
wdfigs       <- paste0(wdproject, "Figures/")

### Chamar funções ----
source(paste0(wdcode, "Assignment scripts/Rfunctions/cosmodist.R")) # Calcula DA_Mpc (Distância angular em Mpc)
source(paste0(wdcode, "Assignment scripts/Rfunctions/poly_fit.R"))

### Decidir se usa distância projetada e distância em z, ou só Yang ----
PM_only <- F # Se F, usa distância projetada e distância em z. Se T, usa só método de Yang.


### Ajustar distância angular para lista de redshifts ----
DA_Mpc_temp <- vector() # Vetor vazio para armazenar as distâncias (angulares?) em Mpc.
zzs         <- seq(0.005, 0.2, 0.005) # Vetor com 40 valores de redshifts, de 0.005 <= z <= 0.2

for(zz in zzs){
  # Usa função cosmodist.R para calcular DA_Mpc:
  DA_Mpc_temp <- append(DA_Mpc_temp, cosmodist(zz, H0 = H_0, Omega_m = Omega_m, Omega_l = Omega_l)$DA_Mpc)
}

DAz_fit <- poly_fit(x = zzs, y = DA_Mpc_temp, degree = 10, nit = 1, nsigma = 50)
plot(zzs, DAz_fit$fitted.values)

source(paste0(wdcode,"Assignment scripts/Rfunctions/NFW.R")) # voltar

### Inicializar data frames de controle info_data ----
info_data             <- data.frame(matrix(ncol = 5, nrow = 0, dimnames = list(NULL, c("order", "info", "nr_gals", "nr_grupos", "nr_centrais"))))
info_data$order       <- as.numeric(info_data$order)
info_data$info        <- as.character(info_data$info)
info_data$nr_gals     <- as.numeric(info_data$nr_gals)
info_data$nr_grupos   <- as.numeric(info_data$nr_grupos)
info_data$nr_centrais <- as.numeric(info_data$nr_centrais)

### Definir tabelas de entrada ----
gals_input_file   <- "SDSS_DR18_Legacy_MGS_QSO.csv"
groups_input_file <- "grupos_lim17.csv"

### Definir tabelas de saída ----
output_file        <- paste0(wdassigndata, "assignment2groups_zmax", zmax, "_Rlim", Rlim, "_Ma", Mhalo_min_group, ".csv")
info_groups_output <- paste0(wdassigndata, "Info/assignment2groups", "_zmax", zmax, "_Rlim", Rlim, "_Ma", Mhalo_min_group, ".infogroups")
info_gals_output   <- paste0(wdassigndata, "Info/assignment2groups", "_zmax", zmax, "_Rlim", Rlim, "_Ma", Mhalo_min_group, ".infogals")

### Ler as tabelas de entrada ----

# Colunas para ler da tabela de galáxias:
gals_colnames <- c('igal', 
                   'ra', 
                   'dec', 
                   'z')

# Colunas para ler da tabela de grupos:
group_colnames <- c('groupID',
                    'centralID',
                    'centralRA',
                    'centralDEC',
                    'centralZ',
                    'Ngals', 
                    'logM100', # logMgroup
                    'logM200', 
                    'R200')

data_gals_i   <- fread(paste0(wddata, "Raw data/SDSS/", gals_input_file), select = gals_colnames)

data_groups_i <- fread(paste0(wddata, "Raw data/Lim+17/", groups_input_file), select = group_colnames)
colnames(data_groups_i)[which(colnames(data_groups_i) == "logM100")]  <- "logMgroup"

### Calcular raios virial (rvir) dos grupos ----
H_z <- sqrt(H_0**2 * (Omega_m * (1 + data_groups_i$centralZ)**3 + 1 - Omega_m)) # Constante de Hubble para um dado z
data_groups_i$rvir <- (10**data_groups_i$logMgroup * 2 * 4.302e-3 * 1e-6 / (100 * H_z**2))**(1/3) # Raio virial do grupo.

# Calcular V_VIR ----
data_groups_i$vvir = sqrt(G * 10**(data_groups_i$logMgroup) * Msun2kg / 
                            (data_groups_i$rvir * Mpc2m)) / 1e3 # Velocidade virial do grupo.

# Seleciona amostra de grupos ----
xx.xx = data_groups_i$Ngals >= Ngals_min_group & # pega os grupos que tem mais ou igual a Ngals_min_group
  data_groups_i$logMgroup >= Mhalo_min_group &   # pega grupos com massa de halo maior ou igual a Mhalo_min_group
  data_groups_i$centralZ >= zmin &              # pega grupos com redshift maior ou igual a zmin
  data_groups_i$centralZ <= zmax                # pega grupos com redshift até zmax

data_groups <- data_groups_i[xx.xx, ] # a partir daqui nunca mais usa o data_groups_i 

groupIDs        <- unique(data_groups$groupID) # IDs dos grupos
groupIDs_select <- groupIDs
Ngroups         <- length(groupIDs) # quantidade de grupos da amostra selecionada

# Bordas do survey ----
cDeltaz <- abs(data_groups$centralZ - zmin) * 299792 / (1 + data_groups$centralZ)
Dzmin  <- cDeltaz * sqrt(DELTA / 2) / data_groups$vvir 

cDeltaz <- abs(data_groups$centralZ - zmax) * 299792 / (1 + data_groups$centralZ)
Dzmax  <- cDeltaz * sqrt(DELTA / 2) / data_groups$vvir

select_groups <- unique(data_groups$groupID[Dzmin >= N_rvir & 
                                              Dzmax >= N_rvir & 
                                              data_groups$groupID %in% groupIDs_select]) # grupos longe das bordas do survey?

groupIDs_select <- groupIDs_select[groupIDs_select %in% select_groups]
Ngroups_select  <- length(groupIDs_select)
print(sprintf('Number of groups (after exluding groups close zmin or zmax): %i', Ngroups_select))

# Seleciona amostra de galáxias (GALAXIES TO BE ASSIGNED) ----
xx.xx = data_gals_i$z >= zmin & # pega galáxias com redshift maior ou igual a zmin
  data_gals_i$z <= zmax # pega galáxias até redshift zmax

data_gals <- data_gals_i[xx.xx, ] # Nunca mais usa data_gals_i
ngals     <- nrow(data_gals)

print(sprintf('Galaxy z range: %6.4f - %6.4f', min(data_gals$z), max(data_gals$z)))

print(sprintf('Group z range (groups far from edges): %6.4f - %6.4f', 
              min(data_groups$centralZ[data_groups$groupID %in% groupIDs_select]), 
              max(data_groups$centralZ[data_groups$groupID %in% groupIDs_select])))

data_groups$DA_Mpc <- predict.lm(DAz_fit, newdata = data.frame(x = data_groups$centralZ)) # distância angular da central
# ANGULAR SIZES OF GROUPS
data_groups$rvir_deg <- data_groups$rvir / data_groups$DA_Mpc * rad2deg # raio virial do grupo em graus?

# Passo 1: Adicionar novas variáveis ----
igal_output       <- vector() #1
ra_output         <- vector() #1
dec_output        <- vector() #2
z_output          <- vector() #3
groupID_output    <- vector() #4
centralID_output  <- vector() #5
centralRA_output  <- vector() #6
centralDEC_output <- vector() #7
centralZ_output   <- vector() #8
Mgroup_output     <- vector() #9
assign_output     <- vector() #10
Rproj_Mpc_output  <- vector() #11
Rproj_rvir_output <- vector() #12
r3D_Mpc_output    <- vector() #13
r3D_rvir_output   <- vector() #14
vlos_vvir_output  <- vector() #15

# LOOP OVER SAMPLE GROUPS ----
print(sprintf('Assigning galaxies to %i groups...........', Ngroups))
progress = round(Ngroups * seq(0, 1, 0.01))
pp = 0

# O RA, DEC e z do grupo é o definido pela central do lim
groupRA  <- data_groups$centralRA
groupDEC <- data_groups$centralDEC
groupZ   <- data_groups$centralZ

# ASSIGNMENT ----
#group <- "224146"

for(group in groupIDs){
  
  i        <- which(data_groups$groupID == group) # pega o id do grupo
  z_i      <- groupZ[i] # pega o redshift do grupo
  z_ik     <- (z_i + data_gals$z) / 2 # médias entre o redshift do grupo e o redshift de todas as galáxias da amostra
  DA_Mpc_i <- data_groups$DA_Mpc[i] # distância angular em Mpc do grupo
  E_z_i    <- sqrt(Omega_m * (1 + z_i)**3 + Omega_l)
  #Rv20     <- N_rvir * data_groups$rvir_deg + N_rvir * data_groups$rvir_deg[i]
  Rv20     <- (N_rvir * data_groups$rvir_deg + N_rvir * data_groups$rvir_deg[i]) + 0.0015 # incerteza de 5 arcsec no céu
  
  # N_rvir Dz_i + N_rvir Dz
  #Dz20 <- N_rvir * (1 + z_i) * (data_groups$vvir + data_groups$vvir[i]) * sqrt(2 / DELTA) / 299792
  Dz20 <- (N_rvir * (1 + z_i) * (data_groups$vvir + data_groups$vvir[i]) * sqrt(2 / DELTA) / 299792) + 0.001 # (aprox. ~ 300 km/s de incerteza)
  # SEARCH RADIUS FOR GALAXIES = N_rvir RVIR_i
  searchRad <- N_rvir * data_groups$rvir_deg[i]
  
  # ANGULAR DISTANCES TO OTHER GROUPS 
  ra_i  <- groupRA[i]
  dec_i <- groupDEC[i]
  
  distGr <- acos(cos((90.0 - dec_i) * deg2rad) * cos((90.0 - groupDEC) * deg2rad) + 
                   sin((90.0 - dec_i) * deg2rad) * sin((90.0 - groupDEC) * deg2rad) * 
                   cos((ra_i - groupRA) * deg2rad)) * rad2deg
  
  distGr[groupRA == ra_i & groupDEC == dec_i] <- 0
  
  # ANGULAR DISTANCES TO GALAXIES
  ra_i  <- groupRA[i]
  dec_i <- groupDEC[i]
  
  distGal <- acos(cos((90.0 - dec_i) * deg2rad) * cos((90.0 - data_gals$dec) * deg2rad) + 
                    sin((90.0 - dec_i) * deg2rad) * sin((90.0 - data_gals$dec) * deg2rad) * 
                    cos((ra_i - data_gals$ra) * deg2rad)) * rad2deg
  
  distGal[data_gals$ra == ra_i & data_gals$dec == dec_i] <- 0
  
  #########################################
  # SINGLE ASSIGNMENT
  #########################################
  # NEARBY GROUPS: 
  #       DISTANCE < (N_rvir RVIR_i + N_rvir RVIR)   
  #  AND  Delta_z  < N_rvir * (v_vir_i + v_vir) * sqrt(2 / DELTA) / c 
  
  nearbyGr <- distGr <= Rv20 & abs(groupZ - z_i) <= Dz20
  NnearbyGr <- length(distGr[nearbyGr])
  
  # DEFINE DATA.FRAME OF NEARBY GROUPS
  data_groups_2NRv <- data_groups[nearbyGr, ]
  
  # NEARBY GALAXIES: DISTANCE < N_rvir RVIR_i
  nearbyGal  <- distGal <= searchRad
  NnearbyGal <- length(distGal[nearbyGal])
  
  if(NnearbyGal == 0){
    write(group, file = paste0(wdassigndata, "Info/NoNearbyGal", "_zmax", zmax, "_Rlim", Rlim, "_Ma", Mhalo_min_group, ".csv"), append = T)
  }else{
    # DEFINE DATA.FRAME FOR NEARBY GALAXIES
    data_gals_NRv = data_gals[nearbyGal, ]
    
    # INITIALIZE VECTORS
    # Passo 2: Adicionando novas variáveis
    
    if(PM_only == T){temp = 0}else{temp = 1e5}
    
    assign_0   <- matrix(temp, ncol = NnearbyGal)
    groupID    <- matrix(-1, ncol = NnearbyGal)
    centralID  <- matrix(-1, ncol = NnearbyGal)
    centralRA  <- matrix(-1, ncol = NnearbyGal)
    centralDEC <- matrix(-1, ncol = NnearbyGal)
    centralZ   <- matrix(-1, ncol = NnearbyGal)
    Mgroup     <- matrix(-1, ncol = NnearbyGal)
    
    vlos_vvir  <- matrix(-1, ncol = NnearbyGal)
    Rproj_Mpc  <- matrix(-1, ncol = NnearbyGal)
    Rproj_rvir <- matrix(-1, ncol = NnearbyGal)
    r3D_Mpc    <- matrix(-1, ncol = NnearbyGal)
    r3D_rvir   <- matrix(-1, ncol = NnearbyGal)
    
    # LOOP OVER NEARBY GROUPS
    # groupRA20Rv  <- data_groups_2NRv$bcgRA
    # groupDEC20Rv <- data_groups_2NRv$bcgDEC
    # groupZ20Rv   <- data_groups_2NRv$bcgZ
    
    # aquiii2
    groupRA20Rv  <- data_groups_2NRv$centralRA
    groupDEC20Rv <- data_groups_2NRv$centralDEC
    groupZ20Rv   <- data_groups_2NRv$centralZ
    
    for(j in 1:NnearbyGr){
      
      z_j      <- groupZ20Rv[j]
      DA_Mpc_j <- data_groups_2NRv$DA_Mpc[j]
      
      # ANGULAR DISTANCES TO NEARBY GALAXIES
      ra_j  <- groupRA20Rv[j]
      dec_j <- groupDEC20Rv[j]
      
      dist_j <- acos(cos((90.0 - dec_j) * deg2rad) * 
                       cos((90.0 - data_gals_NRv$dec) * deg2rad) + 
                       sin((90.0 - dec_j) * deg2rad) * 
                       sin((90.0 - data_gals_NRv$dec) * deg2rad) * 
                       cos((ra_j - data_gals_NRv$ra) * deg2rad)) * rad2deg
      
      dist_j[data_gals_NRv$ra == ra_j & data_gals_NRv$dec == dec_j] <- 0
      
      # ASSIGNMENT EQUATION
      
      # ----------------------
      # R >= Rlim r_vir
      # ----------------------
      R             <- dist_j * deg2rad * DA_Mpc_j / data_groups_2NRv$rvir[j]
      cDeltaz       <- abs(z_j - data_gals_NRv$z) * 299792 / (1 + z_j)
      Dz            <- cDeltaz * sqrt(DELTA / 2) / data_groups_2NRv$vvir[j]
      assign_zSpace <- sqrt(R**2 + Dz**2) 
      
      # ----------------------
      # R < Rlim r_vir
      # ----------------------
      A <- 10**(0.52 + (0.905 - 0.52) * exp(-0.617 * z_j**1.21)) # Dutton & Maccio, 2014, MNRAS, 441, 3359
      B <- -0.101 + 0.026 * z_j                                  # Dutton & Maccio, 2014, MNRAS, 441, 3359
      
      c200        <- A * (10**data_groups_2NRv$logM200[j] * (H_0/100) / 1e12)**B
      conc        <- data_groups_2NRv$rvir[j] / (data_groups_2NRv$R200[j] / c200)
      sigma_group <- data_groups_2NRv$vvir[j] * eta
      
      Density_NFW <- getSigma_NFW(R * data_groups_2NRv$rvir[j], 
                                  data_groups_2NRv$rvir[j], 
                                  conc, 10**data_groups_2NRv$logMgroup[j]) # NFW.R
      
      cDeltaz <- abs(z_j - data_gals_NRv$z) * 299792
      
      pGauss <- (1/sqrt(2*pi)) * 
        (299792 / sigma_group * (1 + z_j)) * 
        exp(-1 * cDeltaz**2 / (2 * sigma_group**2 * (1 + z_j)**2))
      
      assign_PM <- (H_0 / 299792) * (Density_NFW / rho_mean_0) * pGauss
      
      # ----------------------
      # NORMALIZATION AT R = Rlim r_vir
      # ----------------------
      H_z <- sqrt(H_0**2 * (Omega_m * (1 + z_j)**3 + 1 - Omega_m))
      aa2 <- -1 / (DELTA * eta**2)
      f1  <- 1 / 3
      
      if(conc > 1){f1 = (1 - acos(1 / (Rlim * conc)) / 
                           sqrt(abs((conc * Rlim)**2 - 1))) / ((conc * Rlim)**2 - 1)}
      if(conc < 1){f1 = (1 - acosh(1 / (conc * Rlim)) / 
                           sqrt(abs((conc * Rlim)**2 - 1))) / ((conc * Rlim)**2 - 1)}
      
      gc <- 1 / (log(conc + 1) - conc/(conc + 1))
      
      AA <- (2 / 3) * (H_z / H_0) * conc**2 * gc / 
        (Omega_m * eta * (1 + z_j)) * sqrt(DELTA / pi) * f1
      
      aa1 <- log(AA) + Rlim**2 / (100 * eta**2)
      aa  <- c(aa1, 0, aa2)
      
      logPM_Yang <- as.function(polynomial(aa))
      PM_lim <- exp(aa[1])
      
      if(PM_only == F){
        # ----------------------
        # ALL R
        # ----------------------
        # Correção da Marina:
        # O problema é essa condição aqui:  assign_temp[R == 0 | assign_PM >= PM_lim] = 0
        # O grupo que passa primeiro por aqui leva as galáxias (inclusive centrais). 
        # Mudei pra manter o assignment do catálogo original quando o código não sabe o que fazer 
        # (pra assign_PM >= PM_lim, o valor não dá pra ser calculado, a função diverge - depois 
        # te explico!)
        
        # assign_temp = assign_0
        # assign_temp[R <= Rlim & R != 0 & assign_PM < PM_lim] = 
        #   sqrt((log(assign_PM[R <= Rlim & R != 0 & assign_PM < PM_lim]) - aa[1]) / aa[3])
        # assign_temp[R > Rlim] = assign_zSpace[R > Rlim]
        # assign_temp[R == 0 | assign_PM >= PM_lim] = 0
        # assignT = assign_temp < assign_0 & assign_temp <= N_rvir 
        
        assign_temp = assign_0
        assign_temp[R <= Rlim & R != 0 & assign_PM < PM_lim] = 
          sqrt((log(assign_PM[R <= Rlim & R != 0 & assign_PM < PM_lim]) - aa[1]) / aa[3])
        assign_temp[R > Rlim] = assign_zSpace[R > Rlim]
        assign_temp[R == 0] = -1
        assign_temp[assign_PM >= PM_lim] = 0
        assign_temp[assign_PM >= PM_lim & 
                      data_gals_NRv$groupID == data_groups_2NRv$groupID[j]] = -1
        # assign_temp[R == 0 | assign_PM >= PM_lim] = 0
        
        assignT = assign_temp < assign_0 & assign_temp <= N_rvir
        
      }else{
        assign_temp = assign_PM
        assignT = assign_PM > assign_0 & assign_PM >= 0.001
      }
      
      #if(data_groups_2NRv$groupID[j] == 1){col = 'red'}else{col = 'black'}
      #points(R[assignT], assign_PM[assignT], col = col, log = 'xy', xlim = c(0.005, 10))
      #print(data_groups_2NRv$groupID[j])
      
      groupID[assignT]    <- data_groups_2NRv$groupID[j]
      Mgroup[assignT]     <- data_groups_2NRv$logMgroup[j]
      centralID[assignT]  <- data_groups_2NRv$centralID[j]
      centralRA[assignT]  <- data_groups_2NRv$centralRA[j]
      centralDEC[assignT] <- data_groups_2NRv$centralDEC[j]
      centralZ[assignT]   <- data_groups_2NRv$centralZ[j]
      
      Rproj_Mpc[assignT]  <- dist_j[assignT] * deg2rad * DA_Mpc_j
      Rproj_rvir[assignT] <- dist_j[assignT] * deg2rad * DA_Mpc_j / data_groups_2NRv$rvir[j]
      
      vlos_vvir[assignT] <- (abs(z_j - data_gals_NRv$z[assignT]) * 299792 / 
                               (1 + z_j)) / data_groups_2NRv$vvir[j]
      
      assign_0[assignT] <- assign_temp[assignT]
      
    } # END for(j in 1:NnearbyGr)
    
    
    # IF THERE ARE ASSIGNMENTS TO GROUP_i, SAVE IN OUTPUT FILE
    
    members <- groupID == data_groups$groupID[i]
    
    # Passo 3: Adicionar novas variáveis
    igal_output      <- c(igal_output, data_gals_NRv$igal[members])
    ra_output         <- c(ra_output,       data_gals_NRv$ra[members])
    dec_output        <- c(dec_output,      data_gals_NRv$dec[members])
    z_output          <- c(z_output,        data_gals_NRv$z[members])
    
    Rproj_Mpc_output  <- c(Rproj_Mpc_output,  Rproj_Mpc[members])
    Rproj_rvir_output <- c(Rproj_rvir_output, Rproj_rvir[members])
    vlos_vvir_output  <- c(vlos_vvir_output,  vlos_vvir[members])
    
    groupID_output    <- c(groupID_output, groupID[members])
    Mgroup_output     <- c(Mgroup_output,  Mgroup[members])
    assign_output     <- c(assign_output,  assign_0[members])
    
    centralID_output  <- c(centralID_output, centralID[members]) 
    centralRA_output  <- c(centralRA_output, centralRA[members])
    centralDEC_output <- c(centralDEC_output, centralDEC[members])
    centralZ_output   <- c(centralZ_output, centralZ[members])
  }
  
  if(pp %in% progress){
    print(paste(sprintf('%i / %i, %4.0f %s -- ', pp, Ngroups, (pp/Ngroups) * 100, '%'), Sys.time()))
  }
  pp = pp + 1
  
} # END for(group in groups)

print('Writing output files..........')

# Passo 4: Adicionar novas variáveis
output_single = data.frame(igal       = igal_output,
                           ra          = ra_output, 
                           dec         = dec_output, 
                           z           = z_output, 
                           groupID     = groupID_output, 
                           centralID   = centralID_output,  
                           centralRA   = centralRA_output,  
                           centralDEC  = centralDEC_output, 
                           centralZ    = centralZ_output,   
                           logMgroup   = Mgroup_output, 
                           assign      = assign_output, 
                           Rproj_Mpc   = Rproj_Mpc_output, 
                           Rproj_rvir  = Rproj_rvir_output,
                           vlos_vvir   = vlos_vvir_output)

#output_single$flag_good = output_single$ra
output_single$flag_good = 0
output_single$flag_good[output_single$groupID %in% groupIDs_select] = 1

info_groups <- data.frame(info = c("tamanho original", "pós-corte", "longe das bordas"),
                          quantidade = c(nrow(data_groups_i), nrow(data_groups), Ngroups_select))

info_gals <- data.frame(info = c("tamanho original", "pós-corte", "pós-assignment", "flag_good==1"),
                        quantidade = c(nrow(data_gals_i), nrow(data_gals), nrow(output_single), length(which(output_single$flag_good == 1))))

## Salvar tabelas ----

write.table(info_groups, info_groups_output, row.names = F)
write.table(info_gals, info_gals_output, row.names = F)
write.csv(output_single, file = output_file, row.names = F)
