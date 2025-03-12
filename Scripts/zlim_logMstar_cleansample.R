##################################################################################################
# PURPOSE:
#    Compute 95% completeness limit for log M_star as a function of z_max
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

read_data = T
library('Hmisc')

input_samples = c('S_BC03', 'S_V15', 'V1_BC03', 'V2_BC03', 'V1_M05', 'V2_M05')
input_samples = 'MAITE'
Lfrac = 0
compute = T

###############################################
# DEFINITIONS
###############################################
# ------------------------
# INPUT FILE
# ------------------------
input_file <- paste0(wddata, 'letter1_sample.csv')

# ------------------------
# z BINS
# ------------------------
z_min <<- 0.005
if(Lfrac > 0){z_max_seq <<- seq(0.02, 0.12, 0.01)}else{z_max_seq <<- seq(0.01, 0.12, 0.01)}

# ------------------------
# READ INPUT FILE
# ------------------------
if(read_data == T){
  print('Reading data.....')
  data_i = read.csv(input_file, header = T, stringsAsFactors = F)
  ngals = nrow(data)
  print('Done.')
}

data_i <- data_i[-which(data_i$lgm_tot_p50 == -9999.000),]
data_i <- na.omit(data_i)

xx.xx = data_i$z <= max(z_max_seq) #& data_i$Lfrac_petro_r >= Lfrac
data = data_i[xx.xx, ]
z = data$z
m_r = data$magPetro_r

for(input_sample in input_samples){
  # ------------------------
  # OUTPUT FIGURE
  # ------------------------
  if(Lfrac == 0){temp = ''}else{temp = paste('Lfrac', Lfrac, sep = '')}
  
  # FIGURES
  ifelse(!dir.exists(paste0(figswd, "zlim_logMstar")), 
         dir.create(paste0(figswd, "zlim_logMstar")), FALSE)
  
  output_fig1 = paste(figswd, 'zlim_logMstar/logM_zlim_', 
                      input_sample, temp, '.pdf', sep = '')
  output_fig2 = paste(figswd, 'zlim_logMstar/logM_zlim2_', 
                      input_sample, temp, '.pdf', sep = '')
  
  # TABLES
  ifelse(!dir.exists(paste0(wddata, "zlim_logMstar")), 
         dir.create(paste0(wddata, "zlim_logMstar")), FALSE)
  
  output_file1 = paste(wddata, 'zlim_logMstar/logM_zlim_', 
                       input_sample, temp, '.csv', sep = '')
  output_file2 = paste(wddata, 'zlim_logMstar/logM_zlim_', 
                       input_sample, temp, '_coeffs.csv', sep = '')
  
  if(input_sample == 'S_BC03'){M_s = data$logMstar_BC03}
  if(input_sample == 'S_V15'){M_s = data$logMstar_Vaz15}
  if(input_sample == 'S_CB17'){M_s = data$logMstar_CB17}
  if(input_sample == 'V1_BC03'){M_s = data$logMstar_vespa_BC03_dust1}
  if(input_sample == 'V2_BC03'){M_s = data$logMstar_vespa_BC03_dust2}
  if(input_sample == 'V1_M05'){M_s = data$logMstar_vespa_M05_dust1}
  if(input_sample == 'V2_M05'){M_s = data$logMstar_vespa_M05_dust2}
  
  if(input_sample == 'MAITE'){M_s = data$lgm_tot_p50} # aqui
  
  
  if(compute == T){
    # ------------------------
    # 95-PERCENTILES
    # ------------------------
    Nbins = 10
    mass_lim = vector(); a = vector(); b = vector(); c = vector(); N = vector()
    for(i in 1:length(z_max_seq)){
      xx.xx = z >= z_min & z <= z_max_seq[i]
      M_s_i = M_s[xx.xx]
      m_r_i = m_r[xx.xx]
      mm_min = quantile(M_s_i[m_r_i > 17.4 & m_r_i < 17.6], 0.95)
      if(Lfrac == 0){
        mm_max = quantile(M_s_i[m_r_i > 14.5 & m_r_i < 15.5], 0.95)
      }else{
        mm_max = quantile(M_s_i[m_r_i > 15.5 & m_r_i < 16.5], 0.95)
      }
      dd = (mm_max - mm_min) / Nbins
      mass.bin = seq(mm_min, mm_max, dd)
      Mass_med = vector()
      compl95 = matrix(nrow = length(mass.bin)-1)
      for(j in 1:(length(mass.bin)-1)){
        compl95[j] = quantile(m_r_i[M_s_i > mass.bin[j] & M_s_i < mass.bin[j + 1]], 0.95)
        Mass_med[j] = max(M_s_i[M_s_i > mass.bin[j] & M_s_i < mass.bin[j + 1]])
      }
      
      # ------------------------
      # FIT 95-PERCENTILES
      # ------------------------
      fit <- lm(Mass_med ~ poly(compl95, degree = 1, raw = T))
      a[i] = as.numeric(fit$coefficients[1])
      b[i] = as.numeric(fit$coefficients[2])
      #c[i] = as.numeric(fit$coefficients[3])
      c[i] = 0
      N[i] = length(M_s_i)
      
      ###############################################
      # FIGURE 1
      ###############################################
      if(i == 1){
        # ------------------------
        # INITIALIZE FIGURE
        # ------------------------
        pdf(output_fig1, height = 6, width = 6)
        par(mfrow = c(1, 1), mar = c(4, 5, 1, 2), cex.axis = 1.6, cex.lab = 1.6)
        par(font = list(family = 'Times')) 
      }
      
      if(length(M_s_i) > 10000){
        kk.kk = sample(c(1:length(M_s_i)), size = 10000)
      }else{
        kk.kk = c(1:length(M_s_i))
      }
      plot(m_r_i[kk.kk], M_s_i[kk.kk], xlim = c(12, 18), ylim = c(6, 12), 
           cex = 0.2, pch = 19, col = rgb(0, 0, 0, 0.3),
           xlab = expression(paste(italic('m')[italic('r')])), 
           ylab = expression(paste('log ', italic('M')['star'])))
      abline(v = 17.77, lwd = 2, lty = 2, col = 'red')
      
      # ------------------------
      # PLOT 95YH PERCENTILES
      # ------------------------
      points(compl95, Mass_med, pch = 4, lwd = 3, cex = 2, col = 'red')
      # ------------------------
      # PLOT BEST FIT
      # ------------------------
      xx = seq(10, 20); yy = predict.lm(fit, newdata = list(compl95 = xx))
      lines(xx, yy, lty = 1, lwd = 2, col = 'red')
      # ------------------------
      # MASS LIMIT
      # ------------------------
      mass_lim[i] =  predict.lm(fit, newdata = list(compl95 = 17.77))
      abline(h = mass_lim[i], lwd = 2, lty = 2, col = 'red')
      
      # ------------------------
      # LEGEND
      # ------------------------
      legend(12, 7.2, c(paste('z_max = ', z_max_seq[i]), 
                        sprintf('log M_lim = %6.3f', mass_lim[i])), bty = 'n', cex = 1.4)
      
      # ------------------------
      # AXIS
      # ------------------------
      axis(1, seq(11, 19, 0.2), labels = F, tcl = -0.3)
      
      axis(3, seq(11, 19, 1), labels = F, tcl = -0.5)
      axis(3, seq(11, 19, 0.2), labels = F, tcl = -0.3)
      
      axis(2, seq(-30, -10, 0.2), labels = F, tcl = -0.3)
      
      axis(4, seq(-30, -10, 0.2), labels = F, tcl = -0.3)
      axis(4, seq(-30, -10, 1), labels = F, tcl = -0.5)
    }
    dev.off()
  }
  ###############################################
  # SAVE DATA 
  ###############################################
  output_data = data.frame(z_max_seq = z_max_seq, logM = mass_lim, Ngals_bin = N)
  write.csv(output_data, file = output_file1, row.names = F)
  
  if(Lfrac == 0){degree_ = 1}else{degree_ = 3}
  log_z_max_seq = log10(z_max_seq)
  fit_z_max_seq <- lm(log_z_max_seq ~ poly(mass_lim, degree = degree_, raw = T))
  output_coeffs = as.numeric(fit_z_max_seq$coefficients)
  write(c("# Best fit coefficients log(z_max) = a0 + a1 * logM + ... + aN * logM^N"), file = output_file2)
  write(output_coeffs, file = output_file2, append = T)
  
  ###############################################
  # FIGURE 2
  ###############################################
  zs_to_plot <<- c(0.08, 0.06, 0.04)
  # ------------------------------------------------
  # COLORS FOR EACH Z_LIM
  # ------------------------------------------------
  alphaa = 255
  cols <- c(rgb(255, 51, 51, alphaa, maxColorValue = 255), 
             rgb(51, 255, 51, alphaa, maxColorValue = 255), 
             rgb(51, 51, 255, alphaa, maxColorValue = 255))
  cols2 <- c(rgb(204, 0, 0, 255, maxColorValue = 255), 
              rgb(0, 204, 0, 255, maxColorValue = 255), 
              rgb(0, 0, 204, 255, maxColorValue = 255))
  cexx = c(0.3, 0.3, 0.3)
  
  # ------------------------------------------------
  # INITIALIZE FIGURE
  # ------------------------------------------------
  cairo_pdf(output_fig2, height = 6.5, width = 7.5)
  par(mfrow = c(1, 1), mar = c(0, 0, 2, 2), oma = c(4, 5, 0, 0), cex.axis = 1.6, cex.lab = 1.8)
  par(font = list(family = "Times")) 
  
  # ------------------------------------------------
  # EXCLUDE POINTS BEHIND THE SUBPLOT
  # ------------------------------------------------
  sub_x = c(12.3, 14.7); sub_y = c(7.2, 9.8)
  m_r_k = m_r[!(m_r <= sub_y[2] & M_s <= sub_x[2])]
  M_s_k = M_s[!(m_r <= sub_y[2] & M_s <= sub_x[2])]
  
  # ------------------------------------------------
  # PLOT ALL SAMPLE
  # ------------------------------------------------
  set.seed(1)
  #kk.kk = sample(c(1:length(M_s_k)), size = length(M_s_k)/10)
  kk.kk = z > max(zs_to_plot) & z < max(z_max_seq)
  zz.zz = sample(c(1:length(M_s_k[kk.kk])), size = length(M_s_k[kk.kk])/25)
  plot(M_s_k[kk.kk][zz.zz], m_r_k[kk.kk][zz.zz], ylim = c(11.7, 18), xlim = c(6.7, 12.2),  pch = 19, cex = 0.3, 
       xlab = '', 
       ylab = '', 
       col = rgb(0.5, 0.5, 0.5, 1), xaxt = 'n', yaxt = 'n')
  abline(v = 17.77, lwd = 2, lty = 2, col = 'black')
  
  title(ylab = expression(paste(italic('r'))), outer = T)
  title(xlab = expression(paste('log (', italic('m'), ' / M'['\u2299'], ')')), outer = T)
  
  # ------------------------------------------------
  # AXIS
  # ------------------------------------------------
  axis(2, seq(11, 19, 1.0), labels = T, tcl = -0.5)
  axis(2, seq(11, 19, 0.2), labels = F, tcl = -0.3)
  axis(4, seq(11, 19, 1.0), labels = F, tcl = -0.5)
  axis(4, seq(11, 19, 0.2), labels = F, tcl = -0.3)
  
  axis(1, seq(0, 20, 1.0), labels = T, tcl = -0.5)
  axis(1, seq(0, 20, 0.2), labels = F, tcl = -0.3)
  axis(3, seq(0, 20, 1.0), labels = F, tcl = -0.5)
  axis(3, seq(0, 20, 0.2), labels = F, tcl = -0.3)
  
  # ------------------------------------------------
  # PLOT SUBSAMPLES WITHIN DIFFERENT Z_MAX
  # ------------------------------------------------
  for(k in 1:length(zs_to_plot)){
    xx.xx = z >= z_min & z <= zs_to_plot[k] & !(m_r <= sub_x[2] & M_s <= sub_y[2])
    m_r_k = m_r[xx.xx]; M_s_k = M_s[xx.xx]
    
    set.seed(1)
    kk.kk = sample(c(1:length(M_s_k)), round(length(M_s_k)/20))
    points(M_s_k[kk.kk], m_r_k[kk.kk], cex = cexx[k], col = cols[k], pch = 19)
    
    # ------------------------------------------------
    # PLOT BEST FIT TO 95th PERCENTILES
    # ------------------------------------------------
    zz.zz = abs(z_max_seq - zs_to_plot[k]) < 0.0001
    xx = seq(10, 20, 0.01); yy = a[zz.zz] + b[zz.zz] * xx + c[zz.zz] * xx**2
    lines(yy, xx, lty = 1, lwd = 2, col = cols2[k])
  }
  abline(h = 17.77, lwd = 2, lty = 2, col = 'black')
  
  for(k in 1:length(zs_to_plot)){
    zz.zz = abs(z_max_seq - zs_to_plot[k]) < 0.0001
    #abline(h = mass_lim[zz.zz], lwd = 2, lty = 2, col = cols2[k])
    points(mass_lim[zz.zz], 17.77, pch = 21, cex = 3, lwd = 3, col = cols2[k])
  }
  
  # ------------------------------------------------
  # LEGEND
  # ------------------------------------------------
  leg3 = expression(paste(italic('z')['max'], ' = 0.04'))
  leg2 = expression(paste(italic('z')['max'], ' = 0.06'))
  leg1 = expression(paste(italic('z')['max'], ' = 0.08'))
  legend(6.5, 17.7, c(leg1, leg2, leg3), pt.cex = 0.6, cex = 1.4, pch = c(19, 19, 19), 
         col = cols2, bg = 'white', box.col = 'white')
  
  # ------------------------------------------------
  # SUBPLOT WITH LOG M_LIM VS. Z_MAX
  # ------------------------------------------------
  mass_lim <<- mass_lim
  ff <- function(mass_lim, z_max_seq){
    par(cex.axis = 1.1, cex.lab = 1.5, mgp = c(1.7, 0.2, 0))
    par(font = list(family = "Times"))
    
    plot(mass_lim, z_max_seq, pch = 21, cex = 1, lwd = 1.5, xlim = c(8, 11.1), ylim = c(0.005, 0.13),
         xlab = expression(paste('log (', italic('m'), ' / M'['\u2299'], ')')), 
         ylab = expression(paste(italic('z')['max'])), log = 'y',
         xaxt = 'n', yaxt = 'n')
    xx = seq(5, 15, 0.01); yy = predict.lm(fit_z_max_seq, newdata = list(mass_lim = xx))
    lines(xx[xx >= min(mass_lim) & xx <= max(mass_lim)], 
          10**yy[xx >= min(mass_lim) & xx <= max(mass_lim)], lwd = 2, col = 'black')
    lines(xx, 10**yy, lwd = 2, lty = 2, col = 'black')
    
    for(k in 1:length(zs_to_plot)){
      zz.zz = abs(z_max_seq - zs_to_plot[k]) < 0.0001
      lines(c(mass_lim[zz.zz], mass_lim[zz.zz]), c(1e-5, z_max_seq[zz.zz]), col = cols2[k], lty = 5, lwd = 1.5)
      lines(c(1e-5, mass_lim[zz.zz]), c(z_max_seq[zz.zz], z_max_seq[zz.zz]), col = cols2[k], lty = 5, lwd = 1.5)
      
      points(mass_lim[zz.zz], z_max_seq[zz.zz], pch = 19, cex = 1, lwd = 1, col = cols2[k])
      points(mass_lim[zz.zz], z_max_seq[zz.zz], pch = 21, cex = 1, lwd = 1.5)
    }
    
    # -----------------------
    # AXIS
    # -----------------------
    axis(1, seq(5, 19, 1.0), labels = T, tcl = 0.4)
    axis(1, seq(5, 19, 0.2), labels = F, tcl = 0.2)
    axis(3, seq(5, 19, 1.0), labels = F, tcl = 0.4)
    axis(3, seq(5, 19, 0.2), labels = F, tcl = 0.2)
    
    axis(2, c(0.01, 0.1), labels = T, tcl = 0.4)
    axis(2, seq(0.01, 0.1, 0.01), labels = F, tcl = 0.2)
    axis(2, seq(0.001, 0.01, 0.001), labels = F, tcl = 0.2)
    axis(4, c(0.01, 0.1), labels = F, tcl = 0.4)
    axis(4, seq(0.01, 0.1, 0.01), labels = F, tcl = 0.2)
    axis(4, seq(0.001, 0.01, 0.001), labels = F, tcl = 0.2)
    
    # -----------------------
    # LEGEND
    # -----------------------
    leg3 = expression(paste(italic('z')['max'], ' = 0.04'))
    leg2 = expression(paste(italic('z')['max'], ' = 0.06'))
    leg1 = expression(paste(italic('z')['max'], ' = 0.08'))
    legend(8.5, 0.018, c('', '', ''), cex = 1, pch = c(19, 19, 19), col = cols2, bty = 'n')
    legend(8.5, 0.018, c(leg1, leg2, leg3), cex = 1, pch = c(21, 21, 21), pt.lwd = 1.5, bty = 'n')
  }
  subplot( ff(mass_lim, z_max_seq), sub_y, sub_x, vadj = 0, hadj = 0)
  
  dev.off()
  
  
  ###############################################
  pdf('Rplot.pdf', height = 4, width = 5)
  par(mfrow = c(1, 1), mar = c(4, 4, 1, 2), cex.axis = 1.4, cex.lab = 1.4)
  par(font = list(family = 'Times')) 
  fit_z_max_seq <- lm(z_max_seq ~ poly(mass_lim, degree = 4, raw = T))
  plot(mass_lim, z_max_seq, pch = 15, cex = 0.8, lwd = 0.8, xlim = c(8.1, 11.1), ylim = c(0.0, 0.12),
       xlab = expression(paste('log', italic('M')['lim'])), 
       ylab = expression(paste(italic('z')['max'])))
  xx = seq(5, 15, 0.01); yy = predict.lm(fit_z_max_seq, newdata = list(mass_lim = xx))
  lines(xx, yy, lwd = 2, col = 'blue')
  
  #xx = seq(5, 15, 0.01); yy = predict.lm(fit_z_max_seq0, newdata = list(mass_lim = xx))
  #lines(xx, yy, lwd = 2, col = 'red', lty = 2)
  
  legend(8.0, 0.12, c('Best fit ~100 z bins', 'Best fit ~20 z bins'), 
         #legend(8.0, 0.12, c('Best fit z_min = 0.005', 'Best fit z_min = 0.01'), 
         lty = c(1, 2), col = c('blue', 'red'), lwd = c(2, 2), bty = 'n', cex = 1.2)
  dev.off()
}