# List the functions in this file
if(!exists('list_functions_mangle')){
  print('Defining functions from mangle_sdss_masks.R:')
  print('-------------------------------------')
  print(' 1. mangle(data_in, dens_arcmin2, NN_out_min, NN_out_max, mangle_dir)')
  print(' 2. mangle2(data_in, mangle_dir)')
  print('    ****** Suggestion: use mangle2 ******')
  print('')
  list_functions_mangle = F
}

###################################################
# IMPLEMENT A FUNCTION WITH THIS:
#---------------------------------------------------
#mangle_input = data.frame(ra, dec, Rdeg); mangle_output = vector()
#NN_f = 1; kk = 1; NN_sub = 20
#while(NN_f < nrow(mangle_input)){
#  NN_i = (kk - 1) * NN_sub + 1
#  NN_f = NN_i + (NN_sub - 1); if(NN_f > nrow(mangle_input)){NN_f = nrow(mangle_input)}
#  
#  print(paste(sprintf('Running Mangle for groups %i to %i of %i groups', NN_i, NN_f, Ngroups), Sys.time()))
#  
#  mangle_input_i = mangle_input[NN_i:NN_f, ]
#  mangle_output_i = mangle(mangle_input_i, dens_arcmin2 = 10, NN_out_max = 10000, verbose = F)
#  mangle_output = c(mangle_output, mangle_output_i)
#  kk = kk + 1
#}
#temp = data.frame(groupID = groupIDs, f_mangle = mangle_output)
#write.csv(temp, file = paste('../data/mangle/', catalogue, '_groupIDs_mangle_f', mangle_minFrac, '_r', N_rvir, '.csv', sep = ''), 
#          row.names = F)

##################################################################################################
# FUNCTION mangle
##################################################################################################
# PURPOSE:
#    
# CALLING SEQUENCE:
#    mangle(data_in, dens_arcmin2, NN_out_min, NN_out_max)
# INPUTS:
#    data_in        data.frame(ra, dec, Rdeg)
#    dens_arcmin2   Density of points [arcmin^(-2)] (DEFAULT = 10)
#    NN_out_min     Minimum number of points (DEFAULT = 100)
#    NN_out_max     Maximum number of points (DEFAULT = 10000)
# PARAMETERS:
#    
# OUTPUT:
#    ff             vector with fractions (completeness)
# REQUIRED SCRIPTS:
#   
##################################################################################################  

mangle <- function(data_in, dens_arcmin2, NN_out_min, NN_out_max, mangle_dir, verbose){
  # Density of points [arcmin^(-2)]
  if(missing(dens_arcmin2)){dens_arcmin2 = 10}
  
  # Minimum number of points 
  if(missing(NN_out_min)){NN_out_min = 100}
  
  # Maximum number of points 
  if(missing(NN_out_max)){NN_out_max = 10000}
  
  if(missing(verbose)){verbose = T}
  
  if(missing(mangle_dir)){mangle_dir = '../../../../Mangle/mangle2.2/masks/test/'}
  
  # SOME CONSTANTS
  deg2rad = 0.0174532925
  rad2deg = 1 / 0.0174532925

  here = getwd()
  if(verbose == T){print(sprintf('Changing directory to %s', mangle_dir))}
  setwd(mangle_dir)
  ###########################
  # RUN MANGLE
  ###########################
  if(verbose == T){print('Running Mangle.........')}
  id_all = vector()
  
  Nobj = nrow(data_in)
  progress = round(Nobj * seq(0, 1, 0.1))
  
  data_in$id = c(1:Nobj)

  set.seed(1)
  for(i in 1:Nobj){
    id_i = data_in$id[i]
    ra_i = data_in$ra[i]
    dec_i = data_in$dec[i]
    rad_i = data_in$Rdeg[i]
    
    # ********** DOES NOT WORK FOR RA_i DEC_i CLOSE TO THE POLES!!! *********
    op <- options(digits.secs = 6)
    options(op)
    
    area_i = pi * (60 * rad_i)**2
    NN_out_i = round(dens_arcmin2 * area_i)
    if(NN_out_i < NN_out_min){NN_out_i = NN_out_min}
    if(NN_out_i > NN_out_max){NN_out_i = NN_out_max}
    NN = 10 * NN_out_i  #round(4 / pi * NN)
    
    ra_rr = runif(NN, min = (ra_i - 4 * rad_i), max = (ra_i + 4 * rad_i))
    dec_rr = runif(NN, min = (dec_i - 2 * rad_i), max = (dec_i + 2 * rad_i))
    ra_rr[dec_rr > 90] = ra_rr[dec_rr > 90] + 180
    ra_rr[ra_rr > 360] = ra_rr[ra_rr > 360] - 360
    ra_rr[ra_rr < 0] = 360 + ra_rr[ra_rr < 0]
    dec_rr[dec_rr > 90] = 180 - dec_rr[dec_rr > 90]

    dist = acos(cos((90.0 - dec_i) * deg2rad) * cos((90.0 - dec_rr) * deg2rad) + 
                  sin((90.0 - dec_i) * deg2rad) * sin((90.0 - dec_rr) * deg2rad) * 
                  cos((ra_i - ra_rr) * deg2rad)) * rad2deg
    
    if(NN_out_i > length(ra_rr[dist <= rad_i])){NN_out_i = length(ra_rr[dist <= rad_i])}
    ra_out = ra_rr[dist <= rad_i][1:NN_out_i]
    dec_out = dec_rr[dist <= rad_i][1:NN_out_i]

    temp = ra_out; temp[] = id_i
    id_all = c(id_all, temp)
    output = data.frame(ra = ra_out, dec = dec_out, id = temp)

    if(i == 1){
      write.table(output, file = 'mangle_temp.csv', row.names = F, col.names = F)
    }else{
      write.table(output, file = 'temp_coords_i.csv', row.names = F, col.names = F, sep = ' ')
      system(paste('mv mangle_temp.csv temp_coords.csv'))
      system(paste('cat temp_coords.csv temp_coords_i.csv >> mangle_temp.csv'))
    }
    
    if(i %in% progress){
      if(verbose == T){print(paste(sprintf('%i / %i, %4.0f %s -- ', i, Nobj, (i/Nobj) * 100, '%'), Sys.time()))}
    }
  }
  
  system(paste('sh polyid_gals.sh sdss_dr72safe0_res6d.pol ', 
               'mangle_temp.csv mangle_out.dat 0 > lixo.txt', sep = ''))
  
  mangle = count.fields('mangle_out.dat')

  ff = vector()
  for(i in 1:Nobj){
    xx.xx = id_all == i
    ff[i] = length(mangle[xx.xx & mangle == 4]) / length(mangle[xx.xx])
  }
  
  #system('rm mangle_temp.csv temp_coords.csv temp_coords_i.csv')
  system('rm mangle_temp.csv')
  if(verbose == T){print('Done!')}
  
  if(verbose == T){print(sprintf('Changing directory to %s', here))}
  setwd(here)
  return(as.numeric(ff))
}


##################################################################################################
# FUNCTION mangle
##################################################################################################
# PURPOSE:
#    
# CALLING SEQUENCE:
#    mangle(data_in, dens_arcmin2, NN_out_min, NN_out_max)
# INPUTS:
#    data_in        data.frame(ra, dec, Rdeg)
#    dens_arcmin2   Density of points [arcmin^(-2)] (DEFAULT = 10)
#    NN_out_min     Minimum number of points (DEFAULT = 100)
#    NN_out_max     Maximum number of points (DEFAULT = 10000)
# PARAMETERS:
#    
# OUTPUT:
#    ff             vector with fractions (completeness)
# REQUIRED SCRIPTS:
#   
##################################################################################################  

mangle2 <- function(data_in, verbose, mangle_dir){
  if(missing(mangle_dir)){mangle_dir = '../../../../Mangle/mangle2.2/masks/test/'}
  
  if(missing(verbose)){verbose = T}
  
  data_in$id = c(1:nrow(data_in))
  
  # SOME CONSTANTS
  deg2rad = 0.0174532925
  rad2deg = 1 / 0.0174532925
  
  # grid_coords_out.csv --> CREATED BY ALBI/Mangle/mangle2.2/masks/test/generate_coords.R
  if(verbose == T){print(sprintf('Reading file %s.......', paste(mangle_dir, 'grid_coords_out.csv', sep = '')))}
  mangle_output = read.csv(paste(mangle_dir, 'grid_coords_out.csv', sep = ''))
  if(verbose == T){print('Done!')}
  
  ra_rr = mangle_output$ra
  dec_rr = mangle_output$dec
  flag_rr = mangle_output$flag
  ###########################
  # RUN MANGLE
  ###########################
  if(verbose == T){print('Running Mangle.........')}
  ff = vector()
  
  Nobj = nrow(data_in)
  progress = round(Nobj * seq(0, 1, 0.001))
  
  for(i in 1:Nobj){
    id_i = data_in$id[i]
    ra_i = data_in$ra[i]
    dec_i = data_in$dec[i]
    rad_i = data_in$Rdeg[i]
    
    dist = acos(cos((90.0 - dec_i) * deg2rad) * cos((90.0 - dec_rr) * deg2rad) + 
                  sin((90.0 - dec_i) * deg2rad) * sin((90.0 - dec_rr) * deg2rad) * 
                  cos((ra_i - ra_rr) * deg2rad)) * rad2deg
    
    flag_i = flag_rr[dist <= rad_i]
    #plot(ra_rr[dist <= rad_i], dec_rr[dist <= rad_i])
    #points(ra_i, dec_i, pch = 3, cex = 3, lwd = 3)
      
    ff[i] = length(flag_i[flag_i > 0]) / length(flag_i)
    #print(c(length(flag_i[flag_i > 0]), length(flag_i)))
    if(i %in% progress){
      if(verbose == T){print(paste(sprintf('%i / %i, %4.1f %s -- ', i, Nobj, (i/Nobj) * 100, '%'), Sys.time()))}
    }
  }
  
  if(verbose == T){print('Done!')}
  return(as.numeric(ff))
}
