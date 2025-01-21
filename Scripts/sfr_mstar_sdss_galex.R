library(data.table)
library(ggplot2)
library(scales)
library(viridis)
library(cowplot)
library(ggpubr)
library(ggExtra)

source("my_theme.R")

setwd("~/Work/Research/env_quenching/scripts-env-quenching")

datawd     <- "~/Work/Research/env_quenching/env_quenching_data/"
figswd     <- "~/Work/Research/env_quenching/env_quenching_figs/"
input_file <- "DR18_Legacy_MGS_GSWLC-X2_MGS_Lim17_Simard"

df <- fread(paste0(datawd, input_file, ".csv"))

sdss  <- data.frame(logssfr = df$logssfr_sdss, logmstar = df$lgm_tot_p50, survey = "SDSS-MGS")
gswlc <- data.frame(logssfr = df$logssf_gswlc, logmstar = df$logMstar, survey = "GSWLC-X2")

sdss$intercept <- -7.85
sdss$slope     <- -0.3

gswlc$intercept <- -6.55
gswlc$slope     <- -0.45

# Colunas adicionais para plots de conferência
sdss$logssfr_sdss  <- sdss$logssfr
gswlc$logssfr_sdss <- sdss$logssfr
sdss$class_sdss <- sdss$logssfr_sdss > (-0.3*(df$lgm_tot_p50 - 10) - 1.85 - 9) 
gswlc$class_sdss <- sdss$logssfr_sdss > (-0.3*(df$lgm_tot_p50 - 10) - 1.85 - 9) 

# Salvar em 1 tabela só:
df <- rbind(sdss, gswlc)
rm(sdss,gswlc)

df$survey <- factor(df$survey, levels = c("SDSS-MGS", "GSWLC-X2"))

ggplot(df, aes(y = logssfr, x = logmstar)) + 
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon") + 
  scale_fill_viridis(option = "G", direction = -1, name = "Level", limits = c(0,1)) +
  scale_x_continuous(breaks = pretty_breaks(n = 10)) + 
  scale_y_continuous(breaks = pretty_breaks(n = 10), labels = comma) + 
  ylab(ssfr_label) + 
  xlab(mstar_label) + 
  facet_wrap(.~survey) + 

  geom_abline(data = subset(df, df$survey == 'SDSS-MGS'), 
              aes(slope = slope, intercept = intercept), 
              color = 'black', linetype = "longdash") + # Knobel+15
  geom_abline(data = subset(df, df$survey == 'GSWLC-X2'),
              aes(slope = slope, intercept = intercept), 
              color = 'black', linetype = "dotted", linewidth = 1.3) + 
  
  my_theme +
  theme(legend.position = c(0.08, 0.21),
        legend.key.size = unit(0.8, 'cm'),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 22),
        legend.box.background = element_rect(colour = "grey20", linewidth = 1),
        legend.margin = margin(r=10,l=5,t=5,b=15))

ggsave(path = figswd, 
       filename = "ssfr_mstar_sdss_gswlc.pdf",
       device = cairo_pdf, width = width_figs, height = height_figs, units = "in", dpi = 600)

# In `ggplot2`, a 2D density plot represents the density of points in a two-dimensional space. The `level` parameter controls the contour lines that are drawn on the plot. 
# Contour lines represent points of equal density, and `level` determines the number of contour lines to draw. By default, `level = 0.5` and `ggplot2` draws contour lines at 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, and 0.95 quantiles of the density. This means that the contour lines represent areas where the density is 5%, 10%, 25%, 50%, 75%, 90%, and 95% of the maximum density.
# You can customize the levels by specifying the `levels` argument when creating the plot. For example, `levels = c(0.1, 0.25, 0.5, 0.75)` will draw contour lines at 0.1, 0.25, 0.5, and 0.75 quantiles of the density.


ggplot(df, aes(y = logssf_gswlc, x = logMstar, color = logssfr_sdss)) + 
  geom_point(alpha = 0.1) + 
  scale_color_viridis(option = "G", direction = -1) +
  scale_x_continuous(breaks = pretty_breaks(n = 10)) + 
  scale_y_continuous(breaks = pretty_breaks(n = 10), labels = comma) + 
  ylab(ssfr_label) + 
  xlab(mstar_label)


df2 <- df[sample(nrow(df), 10000), ]

ggplot(df, aes(y = logssfr, x = logmstar, color = class_sdss)) + 
  geom_point(alpha = 0.1) + 
  #scale_color_viridis(option = "G", direction = -1) +
  scale_x_continuous(breaks = pretty_breaks(n = 10)) + 
  scale_y_continuous(breaks = pretty_breaks(n = 10), labels = comma) + 
  ylab(ssfr_label) + 
  xlab(mstar_label) + 
  facet_wrap(.~survey) + 
  
  geom_abline(data = subset(df, df$survey == 'SDSS-MGS'), 
              aes(slope = slope, intercept = intercept), 
              color = 'black', linetype = "longdash") + # Knobel+15
  geom_abline(data = subset(df, df$survey == 'GSWLC-X2'),
              aes(slope = slope, intercept = intercept), 
              color = 'black', linetype = "dotted", linewidth = 1.3) + 
  
  my_theme +
  theme(legend.position = c(0.08, 0.21),
        legend.key.size = unit(0.8, 'cm'),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 22),
        legend.box.background = element_rect(colour = "grey20", linewidth = 1),
        legend.margin = margin(r=10,l=5,t=5,b=15))



