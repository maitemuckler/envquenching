catalogue_op = "MPA-JHU" 
logMstar_name = 'logMstar'
logSSFR_name = "logSSFR"

setwd("~/Work/Research/env_quenching/scripts-env-quenching")

read_data = T

library(polynom)
library(bbmle)

#########################################
# OPTIONS
#########################################
# --------------------------
# COMPUTE?
# --------------------------
compute = T


# catalogue_op = "SALIM"
# Rlim            <- 2.5    # assignment = PM Yang for R < Rlim and z-space for R > Rlim
# absMag_lim      <- -12 # -20.4 -18.78
# Mhalo_min_group <- 12.3  # Ma # Group halo masses are greater than Mhalo_min_group
# N_rvir          <- 20 
# z_min  <- 0.01           # z_min, z_max for the galaxy sample
# z_max  <- 0.03
# 
# if(catalogue_op == "SALIM"){logMstar_name = 'logMstar'}
# if(catalogue_op == "SDSS"){logMstar_name = 'lgm_tot_p50'}
# 
# if(catalogue_op == "SALIM"){logSSFR_name = 'logssfr_gswlc'}
# if(catalogue_op == "SDSS"){logSSFR_name = 'specsfr_tot_p50'}

# --------------------------
# GROUP SAMPLE
# --------------------------
catalogue = 'yang'     # yang, gga_yang, maggie, or mock
sample_maggie = 0      # 0 OR 1, used if catalogue = 'gga_yang' OR 'maggie'
Ngals_min_group = 1    # Groups have at least Ngals_min_group galaxies   

compute_Ngals = F  

flag_mangle = F        # Select only groups far from the edges of the Survey?
run_mangle = F         # Run Mangle to select groups far from the edges of the Survey?
# IF flag_mangle = T AND run_mangle = F --> READ PRE-EXISTING FILE
mangle_minFrac = 0.95

z_min = 0.01          # z_min, z_max for the galaxy sample

DELTA = 100
eta = 0.65             # eta = sigma / v_vir

# MOCK
real_space = F         # Used only if catalogue = 'mock'
H_0_mock = 67.3
z_snap_mock = 0.00026
file_bug_mock = 'groups_mock_BUG.csv'

# --------------------------
# ASSIGNMENT SCHEME
# --------------------------
PM_only = F             # Assign galaxies using only Yang's approach (PM_only = T) 
#  OR Yang for R < Rlim r_vir e z-space distance for R > Rlim r_vir (PM_only = F) 
flag_vel = 'vvir'       # sigma_group OR vvir

center = 'L'            # L, M, LW -- WARNING: IF 'M', BE AWARE THAT NOT ALL GALAXIES HAVE THE STELLAR MASS AVAILABLE
#             THOSE WITH NO STELLAR MASSES ARE EXCLUDED IF CENTER = 'M'

#########################################
# INPUT/OUTPUT FILES
#########################################
if(catalogue == 'yang'){
  #input_file = '../data/FinalTable_MasterTable_selectedPars_Vaz15_YangNEW_zoo_extraGals.csv'
  # input_file = paste0(datawd, 'FinalTable_MasterTable_selectedPars_Vaz15_LIM17_zoo_clean-0.1_SALIM_radii.csv')
  input_file = "~/Work/Research/env_quenching/env_quenching_data/DR18_Legacy_MGS_GSWLC-X2_MGS_Simard11_Lim17_clean.csv"
}
if(catalogue == 'gga_yang'){
  input_file = paste('../../../../MAGGIE/Diego/MAGGIE_L_giYangB10_SDSS_BadStripesFalseStarsFalse', 
                     sample_maggie, '_Neff0.0.csv', sep = '')
}
if(catalogue == 'maggie'){
  input_file = paste('../../../../MAGGIE/Diego/MAGGIE_L_SDSS_BadStripesFalseStarsFalse', 
                     sample_maggie, '_Neff0.0.csv', sep = '')
}
if(catalogue == 'mock'){
  input_file = '../../../../MOCKS/Henriques+15/distances_mock_BGGinfo.csv'
}

# input_file_gals = '../data/FinalTable_MasterTable_selectedPars_Vaz15_Yang_zoo_clean.csv'
# input_file_gals = paste0(datawd, 'FinalTable_MasterTable_selectedPars_Vaz15_LIM17_zoo_clean-0.1_SALIM_radii.csv')
input_file =  "~/Work/Research/env_quenching/env_quenching_data/DR18_Legacy_MGS_GSWLC-X2_MGS_Simard11_Lim17_clean.csv"

# OUTPUT FILES
# MAGGIE/GGA samples? 
if(catalogue == 'mock' | catalogue == 'yang'){sample_maggie = ''}  
# mock 3D??? If so, output file name with '_3D'
if(catalogue == 'mock' & real_space == T){temp = '_3D'}else{temp = ''}  
# Mabs limit, z limit, minimum group Mhalo 
temp = paste(temp, '_Mabs', absMag_lim, '_zlim', z_max, '_Ma', Mhalo_min_group, sep = '')
# Ngals limit (appear in the file name if Ngals_min_group > 1)
if(Ngals_min_group > 1){temp = paste(temp, '_Nmin', Ngals_min_group, sep = '')}
# Maximum distance, and flag_vel (sigma_group or v_vir)
temp = paste(temp, '_', N_rvir, 'rvir_', flag_vel, '_center', center, sep = '')

# OUTPUT FILES ----

output_file       <- paste('~/assignment2groups/', Rlim, 'rvir/assignment2groups_', 
                           catalogue, sample_maggie, temp, '_', catalogue_op, '.csv', sep = '')
output_file_multi <- paste('~/assignment2groups/', Rlim, 'rvir/assignment2groups_',
                          catalogue, sample_maggie, temp, '_multi_', catalogue_op, '.csv', sep = '')

# Output directory exists? If not, create it
check_dir <- dir.exists(paste('~/assignment2groups/', Rlim, 'rvir', sep = ''))

if(!check_dir){
  dir.create(paste('~/assignment2groups/', Rlim, 'rvir', sep = ''), recursive = T)
}

print(paste('catalogue: ', catalogue))
print(paste('input file:', input_file))

#########################################
# SOME CONSTANTS
#########################################
deg2rad = 0.0174532925
rad2deg = 1 / 0.0174532925
c = 299792 # km s-1
H_0 = 70.2   # km s-1 Mpc-1
G = 4.302e-9 # Mpc Msun-1 (km s-1)^2
Omega_m = 0.3
Omega_l = 0.7
rho_mean_0 = Omega_m * 3 * H_0**2 / (8 * pi * G) # Msun Mpc-3

G       = 6.67408e-11 # m3 kg-1 s-2
Msun2kg = 1.98855e30 # kg
Mpc2m   = 3.086e22   # m

#########################################
# FUNCTIONS
#########################################
source('Rfunctions/cosmodist.R')
source('Rfunctions/poly_fit.R')
source('Rfunctions/NFW.R')
source('Rfunctions/mangle_SDSS_masks.R')
source('Rfunctions/interp.R')
source('Rfunctions/getR200.R')

# --------------------------
# FUNCTION superpose.eb
# --------------------------
superpose.eb <- function(x, y, ebl, ebu = ebl, length = 0.02, col = col, lwd = 1.5){
  temp = (y - ebl); temp[temp <= 0] = 1e-10
  arrows(x, temp, x, y + ebu, angle = 90, code = 3, length = length, col = col, lwd = lwd)
}

# --------------------------
# FIT DISTANCE
# --------------------------
DA_Mpc_temp = vector()
zzs = seq(0.005, 0.2, 0.005)
i = 1
for(zz in zzs){
  temp <- cosmodist(zz, H0 = H_0, Omega_m = Omega_m, Omega_l = Omega_l)
  DA_Mpc_temp[i] = temp$DA_Mpc
  i = i + 1
}
DAz_fit <- poly_fit(zzs, DA_Mpc_temp, degree = 10, nit = 1, nsigma = 50)

if(compute == T){
  #########################################
  # READ DATA
  #########################################
  if(read_data == T){
    print(paste(Sys.time(), ' - Reading data.....'))
    # --------------------------
    # COLUMNS TO BE READ
    # --------------------------
    if(DELTA == 100){sigma_group_col_name = 'sigma100'}else{sigma_group_col_name = 'sigma200'}
    
    pars_yang = c('ra', 
                  'dec', 
                  'R200_L', 
                  'logM200_L', 
                  'logM100_L', 
                  'Ngals_L', 
                  'z', 
                  'absPetro_r', 
                  'igal_yang', 
                  'groupZ_L', 
                  'groupID_L', 
                  logMstar_name,
                  'groupRA_L', 
                  'groupDEC_L', 
                  logSSFR_name,
                  'velDisp')
    
    pars_yang_new = c('ra', 
                      'dec', 
                      'rvir', 
                      'logMgroup', 
                      'logM100_L', 
                      'Ngals', 
                      'z', 
                      'absPetro_r', 
                      'igal_yang',
                      'groupZ', 
                      'groupID', 
                      'logMstar', 
                      'groupRA', 
                      'groupDEC', 
                      'logSSFR', 
                      'sigma_gal')
    
    pars_maggie = c('id_gal', 
                    'RA_gal', 
                    'DEC_gal', 
                    'redshift_gal', 
                    'stellarmass_gal', 
                    'ssfr_gal', 
                    'is_central_gal',
                    'alpha_group', 
                    'delta_group', 
                    'redshift_group', 
                    'velocity_dispersion_group',
                    'halo_mass_group', 
                    'rvir_group', 
                    'vvir_group', 
                    'flag_good_group',
                    'effective_number_group', 
                    'absmag_r_gal', 
                    'id_group')
    
    pars_maggie_new = c('igal', 
                        'ra', 
                        'dec', 
                        'z', 
                        'logMstar', 
                        'logSSFR', 
                        'is_central_gal', 
                        'groupRA', 
                        'groupDEC', 
                        'groupZ', 
                        #'sigma_group',
                        'logMgroup', 
                        'rvir', 
                        'vvir', 
                        'flag_good_group',
                        'Ngals', 
                        'absPetro_r', 
                        'groupID')
    
    pars_mock = c('ra', 
                  'dec', 
                  'r_crit200', 
                  'm_crit200', 
                  'GroupSize', 
                  'z',
                  'absMag_rDust', 
                  'galaxyID',
                  'ra_bgg', 
                  'dec_bgg', 
                  'z_bgg', 
                  'groupID', 
                  'Mstar', 
                  'logsSFR',
                  'x_gal_Mpc', 
                  'y_gal_Mpc', 
                  'z_gal_Mpc', 
                  'x_group_Mpc', 
                  'y_group_Mpc', 
                  'z_group_Mpc')
    
    pars_mock_new = c('ra', 
                      'dec', 
                      'rvir', 
                      'logMgroup', 
                      'Ngals', 
                      'z', 
                      'absPetro_r', 
                      'igal', 
                      'groupRA', 
                      'groupDEC', 
                      'groupZ', 
                      'groupID',
                      'logMstar', 
                      'logSSFR',
                      'x_gal_Mpc', 
                      'y_gal_Mpc', 
                      'z_gal_Mpc', 
                      'x_group_Mpc', 
                      'y_group_Mpc', 
                      'z_group_Mpc')
    
    pars = c(pars_maggie, pars_yang, pars_mock)
    temp = read.csv(input_file, nrows = 1)
    nn = ncol(temp)
    cols2read = rep("NULL", nn)
    cols2read[names(temp) %in% pars] = NA
    
    # --------------------------
    # READ FILE
    # --------------------------
    data_groups_i = read.csv(input_file, stringsAsFactors = T, colClasses = cols2read)
    if(!('igal' %in% names(data_groups_i))){data_groups_i$igal = c(1:nrow(data_groups_i))}
    
    if(catalogue != 'yang' & catalogue != 'mock'){
      pars = c(pars_yang)
      temp = read.csv(input_file_gals, nrows = 1)
      nn = ncol(temp)
      cols2read = rep("NULL", nn)
      cols2read[names(temp) %in% pars] = NA
      data_gals_i = read.csv(input_file_gals, stringsAsFactors = T, colClasses = cols2read) 
    }
    # --------------------------
    # RENAME COLUMNS
    # --------------------------
    if(catalogue == 'yang'){
      pars = data.frame(pars = pars_yang, pars_new = pars_yang_new, stringsAsFactors = F)
    }
    if(catalogue == 'maggie' | catalogue == 'gga_yang'){
      data_groups_i$alpha_group = data_groups_i$alpha_group * rad2deg
      data_groups_i$delta_group = data_groups_i$delta_group * rad2deg
      data_groups_i$ssfr_gal = data_groups_i$ssfr_gal - 9
      
      pars = data.frame(pars = pars_maggie, pars_new = pars_maggie_new, stringsAsFactors = F)
    }
    if(catalogue == 'mock'){
      pars = data.frame(pars = pars_mock, pars_new = pars_mock_new, stringsAsFactors = F)
      
      data_groups_i$galaxyID = c(1:nrow(data_groups_i))
      #data_groups_i$logSSFR = data_groups_i$galaxyID; data_groups_i$logSSFR[] = NA
    }
    for(pp in 1:nrow(pars)){
      names(data_groups_i)[names(data_groups_i) == pars$pars[pp]] = pars$pars_new[pp]
    }
    if(catalogue == 'yang' | catalogue == 'mock'){
      data_gals_i = data_groups_i
    }else{
      pars = data.frame(pars = pars_yang, pars_new = pars_yang_new, stringsAsFactors = F)
      names(data_gals_i)[names(data_gals_i) == pars$pars[pp]] = pars$pars_new[pp]
    }
    
    # --------------------------
    # COMPUTE R100, M100
    # --------------------------
    if(DELTA == 100){
      data_groups_i$logM200 = data_groups_i$logMgroup
      data_groups_i$r200 = data_groups_i$rvir
      
      if(catalogue == 'yang'){
        data_groups_i$logMgroup = data_groups_i$logM100_L
        H_z = sqrt(H_0**2 * (Omega_m * (1 + data_groups_i$groupZ)**3 + 1 - Omega_m))
      }else{
        print(paste(Sys.time(), ' - Computing M100, R100 from M200, R200.....'))
        M200_grid = seq(min(data_groups_i$logMgroup, na.rm = T), max(data_groups_i$logMgroup, na.rm = T) + 0.1, 0.1)
        
        if(catalogue == 'mock'){
          # Mock snapshot at z = z_snap_mock 
          xx100 = vector(); xx200 = vector()
          
          for(xx in M200_grid){
            temp <- Delta200to100_NFW(xx, z = z_snap_mock)
            xx100 = c(xx100, temp$logM100)
            xx200 = c(xx200, xx)
          }
          
          fit_200to100 <- lm(xx100 ~ poly(xx200, degree = 2, raw = T))
          data_groups_i$logMgroup = predict.lm(fit_200to100, newdata = list(xx200 = data_groups_i$logMgroup))
          
          H_z = sqrt(H_0**2 * (Omega_m * (1 + z_snap_mock)**3 + 1 - Omega_m))
        }else{
          zz_grid = seq(min(data_groups_i$groupZ, na.rm = T), max(data_groups_i$groupZ, na.rm = T), 0.001)
          xx100 = vector(); zz200 = vector(); xx200 = vector()
          
          for(zz in zz_grid){
            for(xx in M200_grid){
              temp <- Delta200to100_NFW(xx, z = zz)
              xx100 = c(xx100, temp$logM100)
              zz200 = c(zz200, zz)
              xx200 = c(xx200, xx)
            }
          }
          
          fit_200to100 <- lm(xx100 ~ poly(xx200, degree = 4, raw = T) + poly(zz200, degree = 2, raw = T))
          data_groups_i$logMgroup = predict.lm(fit_200to100, newdata = list(xx200 = data_groups_i$logMgroup, zz200 = data_groups_i$groupZ))
          
          H_z = sqrt(H_0**2 * (Omega_m * (1 + data_groups_i$groupZ)**3 + 1 - Omega_m))
        }
      }
      # G = 4.302e-3 * 1e-6  -->  Mpc * Msun^-1 * (km/s)^2
      data_groups_i$rvir = (10**data_groups_i$logMgroup * 2 * 4.302e-3 * 1e-6 / (100 * H_z**2))**(1/3) 
    }
    
    # --------------------------
    # COMPUTE V_VIR
    # --------------------------
    data_groups_i$vvir = sqrt(G * 10**(data_groups_i$logMgroup) * Msun2kg / (data_groups_i$rvir * Mpc2m)) / 1e3
    
    # --------------------------
    # COMPUTE N_GALS
    # --------------------------
    if(compute_Ngals == T){
      print(paste(Sys.time(), ' - Computing Ngals.....'))
      groupID = data_groups_i$groupID[data_groups_i$absPetro_r <= absMag_lim]
      temp = as.data.frame(table(groupID)); names(temp)[names(temp) == 'Freq'] = 'Ngals'
      data_groups_i$Ngals_old = data_groups_i$Ngals
      data_groups_i = data_groups_i[, names(data_groups_i) != 'Ngals']
      data_groups_i = merge(data_groups_i, temp, by = "groupID", all.x = T)
      data_groups_i$Ngals[is.na(data_groups_i$Ngals)] = 0
    }
    print(paste(Sys.time(), ' - Done.'))
  }
  
  #########################################
  # SELECT SAMPLE
  #########################################
  
  # SELECT GROUP SAMPLE
  xx.xx = data_groups_i$Ngals >= Ngals_min_group & 
    data_groups_i$logMgroup >= Mhalo_min_group & 
    data_groups_i$groupZ >= z_min & 
    data_groups_i$groupZ <= z_max & 
    !is.na(data_groups_i$groupID)
  
  #if(catalogue == 'mock'){
  #  xx.xx = xx.xx & !(data_groups_i$groupID %in% c(25517, 7763, 23406, 135571, 140983, 152974))  # ???? --> groups with problem
  #}
  if(catalogue == 'gga_yang' | catalogue == 'maggie'){
    xx.xx = xx.xx & data_groups_i$flag_good_group == 'True'
  }
  if(center == 'M'){
    xx.xx = xx.xx & !is.na(data_groups_i$logMstar) & data_groups_i$logMstar > 0
  }
  
  data_groups = data_groups_i[xx.xx, ]
  groupIDs = unique(data_groups$groupID)#; groupIDs = groupIDs[1:10]; data_groups = data_groups[data_groups$groupID %in% groupIDs, ]
  groupIDs_select = groupIDs
  Ngroups = length(groupIDs)
  
  print(sprintf('Number of groups: %i', Ngroups))
  
  # --------------------------
  # SDSS MASKS (MANGLE)
  # --------------------------
  if(flag_mangle == T & run_mangle == T){
    temp = data_groups[, names(data_groups) %in% c('rvir', 'groupID', 'groupRA', 'groupDEC', 'groupZ')]
    temp = unique(temp)
    temp$DA_Mpc = predict.lm(DAz_fit, newdata = data.frame(x = as.numeric(temp$groupZ)))
    Rdeg = (temp$rvir * N_rvir / temp$DA_Mpc) * rad2deg
    ra = temp$groupRA; dec = temp$groupDEC
    
    temp$R20deg = Rdeg
    #write.csv(temp, file = '../data/mangle/sample_mockHenriques_Mmin12_0.01z0.1.csv', row.names = F)
    
    mangle_input = data.frame(ra, dec, Rdeg); mangle_output = vector()
    NN_f = 1; kk = 1; NN_sub = 20
    while(NN_f < nrow(mangle_input)){
      NN_i = (kk - 1) * NN_sub + 1
      NN_f = NN_i + (NN_sub - 1); if(NN_f > nrow(mangle_input)){NN_f = nrow(mangle_input)}
      
      print(paste(sprintf('Running Mangle for groups %i to %i of %i groups', NN_i, NN_f, Ngroups), Sys.time()))
      
      mangle_input_i = mangle_input[NN_i:NN_f, ]
      #mangle_output_i = mangle(mangle_input_i, dens_arcmin2 = 10, NN_out_max = 10000, verbose = F)
      mangle_output_i = mangle2(mangle_input_i, verbose = F)
      mangle_output = c(mangle_output, mangle_output_i)
      kk = kk + 1
    }
    temp = data.frame(groupID = groupIDs, f_mangle = mangle_output)
    write.csv(temp, file = paste(datawd, '/mangle/', catalogue, '_zlim', z_max, '_groupIDs_f_mangle', '_r', N_rvir, '.csv', sep = ''), 
              row.names = F)
  }
  if(flag_mangle == T){
    mangle_output = read.csv(paste(datawd, '/mangle/', catalogue, '_zlim', z_max, '_groupIDs_f_mangle', '_r', N_rvir, '.csv', sep = ''), 
                             stringsAsFactors = F)
    select_groups = mangle_output$groupID[mangle_output$f_mangle >= mangle_minFrac]
    
    #data_groups = data_groups[data_groups$groupID %in% select_groups, ]
    groupIDs_select = groupIDs_select[groupIDs_select %in% select_groups]
    Ngroups_select = length(groupIDs_select)
    
    print(sprintf('Number of groups (after exluding groups close to the edge of the survey): %i', Ngroups_select))
  }
  
  # --------------------------
  # BORDERS OF THE SURVEY
  # --------------------------
  if(catalogue == 'mock' & real_space == T){
    cosmo = cosmodist(z_min, H_0, Omega_m, Omega_l)
    DC_Mpc_min = cosmo$DA_Mpc * (1 + z_min)
    cosmo = cosmodist(z_max, H_0, Omega_m, Omega_l)
    DC_Mpc_max = cosmo$DA_Mpc * (1 + z_max)
    
    dx_i = data_groups$x_group_Mpc * (H_0_mock / H_0)
    dy_i = data_groups$y_group_Mpc * (H_0_mock / H_0)
    dz_i = data_groups$z_group_Mpc * (H_0_mock / H_0)
    distGr_3D_min = abs(sqrt(dx_i**2 + dy_i**2 + dz_i**2) - DC_Mpc_min) / data_groups$rvir
    distGr_3D_max = abs(sqrt(dx_i**2 + dy_i**2 + dz_i**2) - DC_Mpc_max) / data_groups$rvir
    
    dd_min = distGr_3D_min
    dd_max = distGr_3D_max
  }else{
    cDeltaz = abs(data_groups$groupZ - z_min) * 299792 / (1 + data_groups$groupZ)
    Dz_min = cDeltaz * sqrt(DELTA / 2) / data_groups$vvir 
    
    cDeltaz = abs(data_groups$groupZ - z_max) * 299792 / (1 + data_groups$groupZ)
    Dz_max = cDeltaz * sqrt(DELTA / 2) / data_groups$vvir 
    dd_min = Dz_min
    dd_max = Dz_max
  }
  select_groups = unique(data_groups$groupID[dd_min >= N_rvir & dd_max >= N_rvir & data_groups$groupID %in% groupIDs_select])
  #data_groups = data_groups[data_groups$groupID %in% select_groups, ]
  #groupIDs = unique(data_groups$groupID)
  groupIDs_select = groupIDs_select[groupIDs_select %in% select_groups]
  Ngroups_select = length(groupIDs_select)
  print(sprintf('Number of groups (after exluding groups close z_min or z_max): %i', Ngroups_select))
  
  # SELECT GALAXY SAMPLE (GALAXIES TO BE ASSIGNED)
  xx.xx = data_gals_i$z >= z_min & data_gals_i$z <= z_max &
    data_gals_i$absPetro_r <= absMag_lim 
  
  if(center == 'M'){
    xx.xx = xx.xx & !is.na(data_gals_i$logMstar) & data_gals_i$logMstar > 0
  }
  
  data_gals = data_gals_i[xx.xx, ]
  ngals = nrow(data_gals)
  
  print(sprintf('Galaxy z range: %6.4f - %6.4f', min(data_gals$z), max(data_gals$z)))
  print(sprintf('Group z range (groups far from edges): %6.4f - %6.4f', min(data_groups$groupZ[data_groups$groupID %in% groupIDs_select]), 
                max(data_groups$groupZ[data_groups$groupID %in% groupIDs_select])))
  
  ###################################################################################
  # BCGs and MOST MASSIVE GALAXIES COORDS
  bcgRA = vector(); bcgDEC = vector(); bcgZ = vector()
  mmgRA = vector(); mmgDEC = vector(); mmgZ = vector()
  DA_Mpc = vector(); igal = vector()
  rvir_deg = vector()
  x_bgg_Mpc = vector(); y_bgg_Mpc = vector(); z_bgg_Mpc = vector()
  
  print('Getting BCG & MMG coordinates......')
  for(i in 1:Ngroups){
    xx.xx = data_groups$groupID == groupIDs[i]
    mags = data_groups$absPetro_r[xx.xx]
    temp = sort(mags, index.return = T)
    bcgRA[i] = data_groups$ra[xx.xx][temp$ix][1]
    bcgDEC[i] = data_groups$dec[xx.xx][temp$ix][1]
    bcgZ[i] = data_groups$z[xx.xx][temp$ix][1]
    
    if(catalogue == 'mock'){
      x_bgg_Mpc[i] = data_groups$x_gal_Mpc[xx.xx][temp$ix][1]
      y_bgg_Mpc[i] = data_groups$y_gal_Mpc[xx.xx][temp$ix][1]
      z_bgg_Mpc[i] = data_groups$z_gal_Mpc[xx.xx][temp$ix][1]
    }
    
    igal[i] = data_groups$igal[xx.xx][temp$ix][1]
    #cosmo = cosmodist(bcgZ[i], H_0, Omega_m, Omega_l)
    DA_Mpc[i] = predict.lm(DAz_fit, newdata = data.frame(x = bcgZ[i]))
    
    # ANGULAR SIZES OF GROUPS
    rvir_deg[i] = data_groups$rvir[xx.xx][temp$ix][1] / DA_Mpc[i] * rad2deg
    
    Mstar = data_groups$logMstar[xx.xx]
    temp = sort(Mstar, index.return = T, decreasing = T)
    mmgRA[i] = data_groups$ra[xx.xx][temp$ix][1]
    mmgDEC[i] = data_groups$dec[xx.xx][temp$ix][1]
    mmgZ[i] = data_groups$z[xx.xx][temp$ix][1]
  }
  
  if(catalogue == 'mock'){
    temp = data.frame(rvir_deg = rvir_deg, DA_Mpc = DA_Mpc, igal = igal,
                      bcgRA = bcgRA, bcgDEC = bcgDEC, bcgZ = bcgZ, 
                      mmgRA = mmgRA, mmgDEC = mmgDEC, mmgZ = mmgZ,
                      x_bgg_Mpc, y_bgg_Mpc, z_bgg_Mpc)
  }else{
    temp = data.frame(rvir_deg = rvir_deg, DA_Mpc = DA_Mpc, igal = igal,
                      bcgRA = bcgRA, bcgDEC = bcgDEC, bcgZ = bcgZ, 
                      mmgRA = mmgRA, mmgDEC = mmgDEC, mmgZ = mmgZ)
  }
  
  data_groups = merge(data_groups, temp, by = "igal", all.y = T)
  print('Done.')
  
  data_groups = data_groups[, !(names(data_groups) %in% c('igal', 'ra', 'dec', 'z', 'logMstar', 'logSSFR',
                                                          'absPetro_r'))]
  data_groups = unique(data_groups)
  ###################################################################################
  
  # Passo 1: Adicionar novas variáveis
  igal_output = vector()
  ra_output = vector()
  dec_output = vector()
  z_output = vector()
  groupID_output = vector()
  Mgroup_output = vector()
  assign_output = vector()
  logMstar_output = vector()
  logSSFR_output = vector()
  Rproj_Mpc_output = vector()
  Rproj_rvir_output = vector()
  r3D_Mpc_output = vector()
  r3D_rvir_output = vector()
  vlos_vvir_output = vector()
  
  sigma_gal_output    <- vector()


  # LOOP OVER SAMPLE GROUPS
  print(sprintf('Assigning galaxies to %i groups...........', Ngroups))
  progress = round(Ngroups * seq(0, 1, 0.01)); pp = 0
  
  if(center == 'L'){
    groupRA = data_groups$bcgRA
    groupDEC = data_groups$bcgDEC
    groupZ = data_groups$bcgZ
  }
  if(center == 'M'){
    groupRA = data_groups$mmgRA
    groupDEC = data_groups$mmgDEC
    groupZ = data_groups$mmgZ
  }
  if(center == 'LW'){
    groupRA = data_groups$groupRA
    groupDEC = data_groups$groupDEC
    groupZ = data_groups$groupZ
  }
  
  init = T; init_file_bug_mock = T
  for(group in groupIDs){
    i = which(data_groups$groupID == group)[1]
    
    z_i = groupZ[i]; z_ik = (z_i + data_gals$z) / 2
    #temp <- cosmodist(z_i, H0 = H_0, Omega_m = Omega_m, Omega_l = Omega_l)
    DA_Mpc_i = data_groups$DA_Mpc[i]
    
    E_z_i = sqrt(Omega_m * (1 + z_i)**3 + Omega_l)
    # N_rvir RVIR_i + N_rvir RVIR
    Rv20 = N_rvir * data_groups$rvir_deg + N_rvir * data_groups$rvir_deg[i]
    if(catalogue == 'mock' & real_space == T){Rv20 = N_rvir * data_groups$rvir + N_rvir * data_groups$rvir[i]}
    
    # N_rvir Dz_i + N_rvir Dz
    Dz20 = N_rvir * (1 + z_i) * (data_groups$vvir + data_groups$vvir[i]) * sqrt(2 / DELTA) / 299792
    
    # SEARCH RADIUS FOR GALAXIES = N_rvir RVIR_i
    searchRad = N_rvir * data_groups$rvir_deg[i]
    if(catalogue == 'mock' & real_space == T){searchRad = N_rvir * data_groups$rvir[i]}
    
    # ANGULAR DISTANCES TO OTHER GROUPS
    ra_i = groupRA[i]
    dec_i = groupDEC[i]
    distGr = acos(cos((90.0 - dec_i) * deg2rad) * cos((90.0 - groupDEC) * deg2rad) + 
                    sin((90.0 - dec_i) * deg2rad) * sin((90.0 - groupDEC) * deg2rad) * 
                    cos((ra_i - groupRA) * deg2rad)) * rad2deg
    distGr[groupRA == ra_i & groupDEC == dec_i] = 0
    
    # ANGULAR DISTANCES TO GALAXIES
    ra_i = groupRA[i]
    dec_i = groupDEC[i]
    distGal = acos(cos((90.0 - dec_i) * deg2rad) * cos((90.0 - data_gals$dec) * deg2rad) + 
                     sin((90.0 - dec_i) * deg2rad) * sin((90.0 - data_gals$dec) * deg2rad) * 
                     cos((ra_i - data_gals$ra) * deg2rad)) * rad2deg
    distGal[data_gals$ra == ra_i & data_gals$dec == dec_i] = 0
    if(catalogue == 'mock'){      
      dx_i = (data_groups$x_bgg_Mpc - data_groups$x_bgg_Mpc[i]) * (H_0_mock / H_0)
      dy_i = (data_groups$y_bgg_Mpc - data_groups$y_bgg_Mpc[i]) * (H_0_mock / H_0)
      dz_i = (data_groups$z_bgg_Mpc - data_groups$z_bgg_Mpc[i]) * (H_0_mock / H_0)
      distGr_3D = sqrt(dx_i**2 + dy_i**2 + dz_i**2)
      
      dx_i = (data_gals$x_gal_Mpc - data_groups$x_bgg_Mpc[i]) * (H_0_mock / H_0)
      dy_i = (data_gals$y_gal_Mpc - data_groups$y_bgg_Mpc[i]) * (H_0_mock / H_0)
      dz_i = (data_gals$z_gal_Mpc - data_groups$z_bgg_Mpc[i]) * (H_0_mock / H_0)
      distGal_3D = sqrt(dx_i**2 + dy_i**2 + dz_i**2) 
    }
    #########################################
    # MULTIPLE ASSIGNMENT
    #########################################
    if(catalogue == 'mock' & real_space == T){
      # ASSIGNMENT EQUATION
      assign_temp = (distGal_3D / data_groups$rvir[i])  
      assignT = assign_temp <= N_rvir 
    }else{
      # ----------------------
      # R >= Rlim r_vir
      # ----------------------
      R = distGal * deg2rad * DA_Mpc_i / data_groups$rvir[i]
      cDeltaz = abs(z_i - data_gals$z) * 299792 / (1 + z_i) 
      Dz = cDeltaz * sqrt(DELTA / 2) / data_groups$vvir[i]
      assign_zSpace = sqrt(R**2 + Dz**2) 
      
      # ----------------------
      # R < Rlim r_vir
      # ----------------------
      A = 10**(0.52 + (0.905 - 0.52) * exp(-0.617 * z_i**1.21)) # Dutton & Maccio, 2014, MNRAS, 441, 3359
      B = -0.101 + 0.026 * z_i                                  # Dutton & Maccio, 2014, MNRAS, 441, 3359
      if(DELTA == 200){
        c200 = A * (10**data_groups$logMgroup[i] * (H_0/100) / 1e12)**B
        conc = c200
        sigma_group = 397.9 * (10**(data_groups$logMgroup[i] + 0.158) / 1e14)**(0.3214)  # Yang et al. 2007
      }else{
        c200 = A * (10**data_groups$logM200[i] * (H_0/100) / 1e12)**B
        conc = data_groups$rvir[i] / (data_groups$r200[i] / c200)
        sigma_group = data_groups$vvir[i] * eta
      }
      if(flag_vel == 'sigma_group'){
        sigma_group = data_groups$sigma_group[i]; if(sigma_group == 0){sigma_group = data_groups$vvir[i] * eta}
      }
      Density_NFW <- getSigma_NFW(R * data_groups$rvir[i], data_groups$rvir[i], conc, 10**data_groups$logMgroup[i])
      
      cDeltaz = abs(z_i - data_gals$z) * 299792
      pGauss = (1/sqrt(2*pi)) * (299792 / sigma_group * (1 + z_i)) * exp(-1 * cDeltaz**2 / (2 * sigma_group**2 * (1 + z_i)**2))
      
      assign_PM = (H_0 / 299792) * (Density_NFW / rho_mean_0) * pGauss
      
      # ----------------------
      # NORMALIZATION AT R = Rlim r_vir
      # ----------------------
      H_z = sqrt(H_0**2 * (Omega_m * (1 + z_i)**3 + 1 - Omega_m))
      aa2 = -1 / (DELTA * eta**2)
      f1 = 1 / 3
      if(conc > 1){f1 = (1 - acos(1 / (Rlim * conc)) / sqrt(abs((conc * Rlim)**2 - 1))) / ((conc * Rlim)**2 - 1)}
      if(conc < 1){f1 = (1 - acosh(1 / (conc * Rlim)) / sqrt(abs((conc * Rlim)**2 - 1))) / ((conc * Rlim)**2 - 1)}
      
      gc = 1 / (log(conc + 1) - conc/(conc + 1))
      AA = (2 / 3) * (H_z / H_0) * conc**2 * gc / (Omega_m * eta * (1 + z_i)) * sqrt(DELTA / pi) * f1
      aa1 = log(AA) + Rlim**2 / (100 * eta**2)
      
      aa = c(aa1, 0, aa2)
      logPM_Yang <- as.function(polynomial(aa))
      PM_lim = exp(aa[1])
      if(PM_only == F){
        # ----------------------
        # ALL R
        # ----------------------
        assign_temp = assign_zSpace; assign_temp[] = -99
        assign_temp[R <= Rlim & R != 0 & assign_PM < PM_lim] = 
          sqrt((log(assign_PM[R <= Rlim & R != 0 & assign_PM < PM_lim]) - aa[1]) / aa[3])
        assign_temp[R > Rlim] = assign_zSpace[R > Rlim]
        assign_temp[R == 0 | assign_PM >= PM_lim] = 0
        
        assignT = assign_temp <= N_rvir 
      }else{
        assignT = assign_PM >= 0.001 & R < N_rvir # exp(logPM_Yang(N_rvir))
      }
    }
    
    Rproj_Mpc_multi = distGal[assignT] * deg2rad * DA_Mpc_i
    Rproj_rvir_multi = distGal[assignT] * deg2rad * DA_Mpc_i / data_groups$rvir[i]
    vlos_vvir_multi = (abs(z_i - data_gals$z[assignT]) * 299792 / (1 + z_i)) / data_groups$vvir[i]
    groupID_multi = Rproj_Mpc_multi; groupID_multi[] = data_groups$groupID[i]
    logMgroup_multi = Rproj_Mpc_multi; logMgroup_multi[] = data_groups$logMgroup[i]
    
    if(init == T){
      output_multi = data.frame(igal = data_gals$igal[assignT], 
                                ra = data_gals$ra[assignT], 
                                dec = data_gals$dec[assignT], 
                                z = data_gals$z[assignT], 
                                groupID = groupID_multi, 
                                logMgroup = logMgroup_multi,
                                assign = assign_temp[assignT], 
                                logMstar = data_gals$logMstar[assignT], 
                                logSSFR = data_gals$logSSFR[assignT],
                                Rproj_Mpc = Rproj_Mpc_multi, 
                                Rproj_rvir = Rproj_rvir_multi,
                                vlos_vvir = vlos_vvir_multi)
      
      if(catalogue == 'mock'){
        output_multi$r3D_Mpc = distGal_3D[assignT]
        output_multi$r3D_rvir = distGal_3D[assignT] / data_groups$rvir[i]
      }
      init = F
    }else{
      temp = data.frame(igal = data_gals$igal[assignT], 
                        ra = data_gals$ra[assignT], 
                        dec = data_gals$dec[assignT], 
                        z = data_gals$z[assignT], 
                        groupID = groupID_multi, 
                        logMgroup = logMgroup_multi,
                        assign = assign_temp[assignT], 
                        logMstar = data_gals$logMstar[assignT], 
                        logSSFR = data_gals$logSSFR[assignT],
                        Rproj_Mpc = Rproj_Mpc_multi, 
                        Rproj_rvir = Rproj_rvir_multi,
                        vlos_vvir = vlos_vvir_multi)
      
      if(catalogue == 'mock'){
        temp$r3D_Mpc = distGal_3D[assignT]
        temp$r3D_rvir = distGal_3D[assignT] / data_groups$rvir[i]
      }
      output_multi = rbind(output_multi, temp)
    }
    
    #########################################
    # SINGLE ASSIGNMENT
    #########################################
    # NEARBY GROUPS: 
    #       DISTANCE < (N_rvir RVIR_i + N_rvir RVIR)   
    #  AND  Delta_z  < N_rvir * (v_vir_i + v_vir) * sqrt(2 / DELTA) / c 
    if(catalogue == 'mock' & real_space == T){
      nearbyGr <- distGr_3D <= Rv20
    }else{
      nearbyGr <- distGr <= Rv20 & abs(groupZ - z_i) <= Dz20
    }
    NnearbyGr = length(distGr[nearbyGr])
    # DEFINE DATA.FRAME OF NEARBY GROUPS
    data_groups_2NRv = data_groups[nearbyGr, ]
    
    # NEARBY GALAXIES: DISTANCE < N_rvir RVIR_i
    if(catalogue == 'mock' & real_space == T){
      nearbyGal <- distGal_3D <= searchRad
    }else{
      nearbyGal <- distGal <= searchRad
    }
    NnearbyGal = length(distGal[nearbyGal])
    if(NnearbyGal == 0){
      if(init_file_bug_mock == T){
        write.csv(data.frame(groupID = group), file = file_bug_mock, row.names = F)
        init_file_bug_mock = F
      }else{
        write(group, file = file_bug_mock, append = T)
      }
    }else{
      # DEFINE DATA.FRAME FOR NEARBY GALAXIES
      data_gals_NRv = data_gals[nearbyGal, ]
      
      # INITIALIZE VECTORS
      # Passo 2: Adicionando novas variáveis
      if(PM_only == T){temp = 0}else{temp = 1e5}
      assign_0 = matrix(temp, ncol = NnearbyGal)
      groupID = matrix(-1, ncol = NnearbyGal)
      Mgroup = matrix(-1, ncol = NnearbyGal)
      vlos_vvir = matrix(-1, ncol = NnearbyGal)
      Rproj_Mpc = matrix(-1, ncol = NnearbyGal)
      igal = matrix(-1, ncol = NnearbyGal)
      Rproj_rvir = matrix(-1, ncol = NnearbyGal)
      r3D_Mpc = matrix(-1, ncol = NnearbyGal)
      r3D_rvir = matrix(-1, ncol = NnearbyGal)

      sigma_gal    <- matrix(-1, ncol = NnearbyGal)
      
      # LOOP OVER NEARBY GROUPS
      if(center == 'L'){
        groupRA20Rv = data_groups_2NRv$bcgRA
        groupDEC20Rv = data_groups_2NRv$bcgDEC
        groupZ20Rv = data_groups_2NRv$bcgZ
      }
      if(center == 'M'){
        groupRA20Rv = data_groups_2NRv$mmgRA
        groupDEC20Rv = data_groups_2NRv$mmgDEC
        groupZ20Rv = data_groups_2NRv$mmgZ
      }
      if(center == 'LW'){
        groupRA20Rv = data_groups_2NRv$groupRA
        groupDEC20Rv = data_groups_2NRv$groupDEC
        groupZ20Rv = data_groups_2NRv$groupZ
      }
      
      for(j in 1:NnearbyGr){
        z_j = groupZ20Rv[j]
        #temp <- cosmodist(z_j, H0 = H_0, Omega_m = Omega_m, Omega_l = Omega_l)
        DA_Mpc_j = data_groups_2NRv$DA_Mpc[j]
        
        # ANGULAR DISTANCES TO NEARBY GALAXIES
        ra_j = groupRA20Rv[j]
        dec_j = groupDEC20Rv[j]
        dist_j = acos(cos((90.0 - dec_j) * deg2rad) * cos((90.0 - data_gals_NRv$dec) * deg2rad) + 
                        sin((90.0 - dec_j) * deg2rad) * sin((90.0 - data_gals_NRv$dec) * deg2rad) * 
                        cos((ra_j - data_gals_NRv$ra) * deg2rad)) * rad2deg
        dist_j[data_gals_NRv$ra == ra_j & data_gals_NRv$dec == dec_j] = 0
        
        if(catalogue == 'mock'){
          dx_j = (data_gals_NRv$x_gal_Mpc - data_groups_2NRv$x_bgg_Mpc[j]) * (H_0_mock / H_0)
          dy_j = (data_gals_NRv$y_gal_Mpc - data_groups_2NRv$y_bgg_Mpc[j]) * (H_0_mock / H_0)
          dz_j = (data_gals_NRv$z_gal_Mpc - data_groups_2NRv$z_bgg_Mpc[j]) * (H_0_mock / H_0)
          distGal_3D = sqrt(dx_j**2 + dy_j**2 + dz_j**2)
        }
        
        # ASSIGNMENT EQUATION
        if(catalogue == 'mock' & real_space == T){
          R = distGal_3D
          assign_temp = distGal_3D / data_groups_2NRv$rvir[j]
          assignT = assign_temp < assign_0 & assign_temp <= N_rvir  
        }else{
          # ----------------------
          # R >= Rlim r_vir
          # ----------------------
          R = dist_j * deg2rad * DA_Mpc_j / data_groups_2NRv$rvir[j]
          cDeltaz = abs(z_j - data_gals_NRv$z) * 299792 / (1 + z_j)
          Dz = cDeltaz * sqrt(DELTA / 2) / data_groups_2NRv$vvir[j]
          assign_zSpace = sqrt(R**2 + Dz**2) 
          
          # ----------------------
          # R < Rlim r_vir
          # ----------------------
          A = 10**(0.52 + (0.905 - 0.52) * exp(-0.617 * z_j**1.21)) # Dutton & Maccio, 2014, MNRAS, 441, 3359
          B = -0.101 + 0.026 * z_j                                  # Dutton & Maccio, 2014, MNRAS, 441, 3359
          if(DELTA == 200){
            c200 = A * (10**data_groups_2NRv$logMgroup[j] * (H_0/100) / 1e12)**B
            conc = c200
            sigma_group = 397.9 * (10**(data_groups_2NRv$logMgroup[j] + 0.158) / 1e14)**(0.3214)  # Yang et al. 2007
          }else{
            c200 = A * (10**data_groups_2NRv$logM200[j] * (H_0/100) / 1e12)**B
            conc = data_groups_2NRv$rvir[j] / (data_groups_2NRv$r200[j] / c200)
            sigma_group = data_groups_2NRv$vvir[j] * eta
          }
          if(flag_vel == 'sigma_group'){
            sigma_group = data_groups_2NRv$sigma_group[j]; if(sigma_group == 0){sigma_group = data_groups_2NRv$vvir[j] * eta}
          }
          Density_NFW <- getSigma_NFW(R * data_groups_2NRv$rvir[j], data_groups_2NRv$rvir[j], 
                                      conc, 10**data_groups_2NRv$logMgroup[j])
          
          cDeltaz = abs(z_j - data_gals_NRv$z) * 299792
          pGauss = (1/sqrt(2*pi)) * (299792 / sigma_group * (1 + z_j)) * exp(-1 * cDeltaz**2 / (2 * sigma_group**2 * (1 + z_j)**2))
          
          assign_PM = (H_0 / 299792) * (Density_NFW / rho_mean_0) * pGauss
          
          # ----------------------
          # NORMALIZATION AT R = Rlim r_vir
          # ----------------------
          H_z = sqrt(H_0**2 * (Omega_m * (1 + z_j)**3 + 1 - Omega_m))
          aa2 = -1 / (DELTA * eta**2)
          f1 = 1 / 3
          if(conc > 1){f1 = (1 - acos(1 / (Rlim * conc)) / sqrt(abs((conc * Rlim)**2 - 1))) / ((conc * Rlim)**2 - 1)}
          if(conc < 1){f1 = (1 - acosh(1 / (conc * Rlim)) / sqrt(abs((conc * Rlim)**2 - 1))) / ((conc * Rlim)**2 - 1)}
          
          gc = 1 / (log(conc + 1) - conc/(conc + 1))
          AA = (2 / 3) * (H_z / H_0) * conc**2 * gc / (Omega_m * eta * (1 + z_j)) * sqrt(DELTA / pi) * f1
          aa1 = log(AA) + Rlim**2 / (100 * eta**2)
          
          aa = c(aa1, 0, aa2)
          logPM_Yang <- as.function(polynomial(aa))
          PM_lim = exp(aa[1])
          
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
            assign_temp[assign_PM >= PM_lim & data_gals_NRv$groupID == data_groups_2NRv$groupID[j]] = -1
            # assign_temp[R == 0 | assign_PM >= PM_lim] = 0
            
            assignT = assign_temp < assign_0 & assign_temp <= N_rvir
            
          }else{
            assign_temp = assign_PM
            assignT = assign_PM > assign_0 & assign_PM >= 0.001
          }
        }
        #if(data_groups_2NRv$groupID[j] == 1){col = 'red'}else{col = 'black'}
        #points(R[assignT], assign_PM[assignT], col = col, log = 'xy', xlim = c(0.005, 10))
        #print(data_groups_2NRv$groupID[j])
        
        groupID[assignT] = data_groups_2NRv$groupID[j]
        Mgroup[assignT] = data_groups_2NRv$logMgroup[j]
        
        Rproj_Mpc[assignT] = dist_j[assignT] * deg2rad * DA_Mpc_j
        Rproj_rvir[assignT] = dist_j[assignT] * deg2rad * DA_Mpc_j / data_groups_2NRv$rvir[j]
        
        vlos_vvir[assignT] = (abs(z_j - data_gals_NRv$z[assignT]) * 299792 / (1 + z_j)) / data_groups_2NRv$vvir[j]
        
        igal[assignT] = data_gals_NRv$igal[assignT]
        assign_0[assignT] = assign_temp[assignT]
        if(catalogue == 'mock'){
          r3D_Mpc[assignT] = distGal_3D[assignT]
          r3D_rvir[assignT] = distGal_3D[assignT] / data_groups_2NRv$rvir[j]
        }
      } # END for(j in 1:NnearbyGr)
      
      # IF THERE ARE ASSIGNMENTS TO GROUP_i, SAVE IN OUTPUT FILE
      
      members <- groupID == data_groups$groupID[i]
      
      # Passo 3: Adicionar novas variáveis
      igal_output = c(igal_output, data_gals_NRv$igal[members])
      ra_output = c(ra_output, data_gals_NRv$ra[members])
      dec_output = c(dec_output, data_gals_NRv$dec[members])
      z_output = c(z_output, data_gals_NRv$z[members])
      
      groupID_output = c(groupID_output, groupID[members])
      Mgroup_output = c(Mgroup_output, Mgroup[members])
      assign_output = c(assign_output, assign_0[members])
      
      logMstar_output = c(logMstar_output, data_gals_NRv$logMstar[members])
      logSSFR_output = c(logSSFR_output, data_gals_NRv$logSSFR[members])
      
      Rproj_Mpc_output = c(Rproj_Mpc_output, Rproj_Mpc[members])
      Rproj_rvir_output = c(Rproj_rvir_output, Rproj_rvir[members])
      vlos_vvir_output = c(vlos_vvir_output, vlos_vvir[members])
      
      sigma_gal_output    <- c(sigma_gal_output, data_gals_NRv$sigma_gal[members])
      
      if(catalogue == 'mock'){
        r3D_Mpc_output = c(r3D_Mpc_output, r3D_Mpc[members])
        r3D_rvir_output = c(r3D_rvir_output, r3D_rvir[members])
      }
    }
    
    if(pp %in% progress){
      print(paste(sprintf('%i / %i, %4.0f %s -- ', pp, Ngroups, (pp/Ngroups) * 100, '%'), Sys.time()))
    }
    pp = pp + 1
    
  } # END for(group in groups)
  
  print('Writing output files..........')
  # Passo 4: Adicionar novas variáveis
  output_single = data.frame(igal = igal_output, 
                             ra = ra_output, 
                             dec = dec_output, 
                             z = z_output, 
                             groupID = groupID_output, 
                             logMgroup = Mgroup_output, 
                             assign = assign_output, 
                             logMstar = logMstar_output, 
                             logSSFR = logSSFR_output,
                             Rproj_Mpc = Rproj_Mpc_output, 
                             Rproj_rvir = Rproj_rvir_output,
                             vlos_vvir = vlos_vvir_output,
                             
                             sigma_gal = sigma_gal_output)
  
  if(catalogue == 'mock'){
    output_single$r3D_Mpc = r3D_Mpc_output
    output_single$r3D_rvir = r3D_rvir_output
  }
  output_single$flag_good = output_single$ra; output_single$flag_good[] = 0
  output_single$flag_good[output_single$groupID %in% groupIDs_select] = 1
  write.csv(output_single, file = output_file, row.names = F)
  
  output_multi$flag_good = output_multi$ra; output_multi$flag_good[] = 0
  output_multi$flag_good[output_multi$groupID %in% groupIDs_select] = 1
  write.csv(output_multi, file = output_file_multi, row.names = F)
  print('Done!')
}
