library(ggplot2)
library(ggh4x)
library(readr)
library(dplyr)

setwd("~/Work/Research/env_quenching")

source("env_quenching_scripts/Themes/my_theme.R")

datawd  <- "env_quenching_data/"
figswd  <- "env_quenching_figs/"

cat_vec         <- c("GSWLC", "MPA-JHU")
zmax_vec        <- 0.1
Mhalo_min_group <- 13
lim_vec         <- c("_nocentral", "")

for (lim_op in lim_vec) {
  for (catalogue_op in cat_vec) {
    for (z_max in zmax_vec) {
      
      median_bins <- read_csv(paste0(datawd,
                                     "output_model/median_bins/",
                                     "median_bins_",
                                     catalogue_op,
                                     "_zmax", z_max,
                                     "_Ma", Mhalo_min_group,
                                     "_flag_good==1_MANGLE",
                                     lim_op,
                                     ".csv")) %>% mutate(Ma = "13.0")
      median_bins$label <- paste0(median_bins$label,"-Ma=" ,median_bins$Ma)
      
      model13 <- read_csv(paste0(datawd,
                                 "output_model/model/",
                                 "model_",
                                 catalogue_op,
                                 "_zmax", z_max,
                                 "_Ma", Mhalo_min_group,
                                 "_flag_good==1_MANGLE",                                
                                 lim_op,
                                 ".csv")) %>% mutate(Ma = "13.0")
      
      model12.3 <- read_csv(paste0(datawd,
                                   "output_model/model/",
                                   "model_",
                                   catalogue_op,
                                   "_zmax", z_max,
                                   "_Ma", "12.3",
                                   "_flag_good==1_MANGLE",                                
                                   lim_op,
                                   ".csv")) %>% mutate(Ma = "12.3")
      
      model12.3 <- subset(model12.3, model12.3$logMgroup != "log[10]~M['h']~(12.3-13.0)")
      model <- rbind(model12.3, model13)
      model$label <- paste0(model$label,"-Ma=" ,model$Ma)
      
      if(z_max == 0.03){
        width_figs  <- 9
        height_figs <- 6
        size_text_facet <- 13
      }
      
      if(z_max == 0.1){
        width_figs  <- 9
        height_figs <- 7
        size_text_facet <- 12
      }
      
      model$AGN_whan <- factor(model$AGN_whan, levels = c("AGN", "Non-AGN"))
      median_bins$AGN_whan <- factor(median_bins$AGN_whan, levels = c("AGN", "Non-AGN"))
      
      p <- ggplot() + 
        geom_line(data = model, aes(x = 10^logRproj, y = pred, 
                                    color = label, linetype = AGN_whan), linewidth = 1) + 
        geom_point(data = median_bins, aes(x = temp_rbins, y = fsf, color = label, shape = AGN_whan), size = 2) +
        
        geom_errorbar(data = median_bins, aes(ymin = confint_lw, ymax = confint_up, 
                                              x = temp_rbins, color = label),
                      width = 0.05, alpha = 0.4) + 
        
        geom_ribbon(data = model, aes(ymin = lwr, ymax = upr, x = 10^logRproj, 
                                      fill = label),
                    inherit.aes = FALSE, alpha = 0.3) +
        
        # Centrais
        geom_point(data = subset(median_bins, median_bins$AGN_whan == "AGN"), 
                   aes(x = 10^-1.15, y = fsf_c, color = label), size = 3, shape = 1) +
        geom_errorbar(data = subset(median_bins, median_bins$AGN_whan == "AGN"), 
                      aes(ymin = confint_lw_c, ymax = confint_up_c, x = 10^-1.15, color = label),
                      width = 0.05, alpha = 0.3) +
        
        geom_point(data = subset(median_bins, median_bins$AGN_whan == "Non-AGN"), 
                   aes(x = 10^-1.26, y = fsf_c, color = label), size = 3, shape = 0) +
        geom_errorbar(data = subset(median_bins, median_bins$AGN_whan == "Non-AGN"), 
                      aes(ymin = confint_lw_c, ymax = confint_up_c, x = 10^-1.26, color = label),
                      width = 0.05, alpha = 0.3) +
        
        facet_nested(logMgroup  ~ logvelDisp_e, labeller = label_parsed) + 
        
        scale_x_log10(n.breaks = 5) + 
        annotation_logticks(sides = "b") + 
        scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
        coord_cartesian(xlim = 10^c(-1.28, 1.35), ylim = c(0,1)) + 
        
        
        scale_color_manual(values = c("#78B3F7", "#052d5a", "#78B3F7", "#052d5a", 
                                      "#FF9999", "#d00000", "#FF9999", "#d00000")) +
        scale_fill_manual(values = c("#78B3F7", "#052d5a", "#78B3F7", "#052d5a", 
                                     "#FF9999", "#d00000", "#FF9999", "#d00000")) +
        scale_linetype_manual(values = c("solid", "dotdash")) + 
        scale_shape_manual(values = c(19, 15)) + 
        
        labs(x = expression(R/r["vir"]), 
             y = expression(f["SF"]), 
             color = expression(log['10']~M['★']),
             fill = expression(log['10']~M['★']),
             shape = "AGN",
             linetype = "AGN") + 
        my_theme + 
        theme(axis.text = element_text(size = 16),
              axis.title = element_text(size = 18),
              legend.text = element_text(size = 14),
              legend.title = element_text(size = 16),
              legend.position = "bottom",
              strip.text = element_text(size = size_text_facet)) + 
        
        guides(fill = guide_legend(nrow = 2, byrow = F),
               color = guide_legend(nrow = 2, byrow = F),
               linetype = guide_legend(nrow = 2, byrow = F),
               shape = guide_legend(nrow = 2, byrow = F),
               legend.spacing.y = unit(1.0, 'cm'))
      
      p 
      
      ggsave(path = paste0(figswd, "/Models/"),
             filename = paste0("logistica_Ma", Mhalo_min_group,
                               "_zmax",z_max,
                               "_",catalogue_op,
                               lim_op,
                               ".pdf"),
             device = cairo_pdf, width = width_figs, height = height_figs, 
             units = "in", dpi = 600)
    }
  }
}

