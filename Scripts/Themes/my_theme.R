source("~/Work/Research/Astronomy/Projects/envquenching/Scripts/Themes/ggplot_theme_Publication-2.R")

my_theme <- theme(panel.grid.major = element_line(colour = "grey90"),
                  panel.grid.minor = element_line(colour = "grey90"),
                  panel.background = element_blank(), 
                  panel.border = element_rect(linetype = "solid", fill = NA),
                  panel.spacing = unit(2, "lines"),
                  
                  axis.line = element_blank(),
                  axis.text = element_text(size = 30, color = "black"),
                  axis.title = element_text(size = 36, color = "black"),
                  
                  plot.title = element_text(size = 36, hjust = 0.5), 
                  plot.tag.position = c(0.05, 0.95),
                  plot.margin = unit(c(1.5,1,1,1), "cm"),
                  
                  plot.subtitle = element_text(size = 28, hjust = 0.5), 
                  
                  legend.position = "bottom",
                  legend.text = element_text(size = 32),
                  legend.title = element_text(size = 36),
                  
                  strip.background = element_rect(colour = "black", fill = "white"),
                  strip.placement = "outside",
                  strip.text = element_text(size = 34),

                  text = element_text(family = "serif"))

my_guides <- guides(fill = guide_legend(nrow = 2, byrow = F),
                    color = guide_legend(nrow = 2, byrow = F),
                    linetype = guide_legend(nrow = 2, byrow = F),
                    shape = guide_legend(nrow = 2, byrow = F),
                    legend.spacing.y = unit(1.0, 'cm'))

# colors:
cores_misclass <- c(TN = "red3", TP = "#023e8a", FN = "#ff8fab", FP = "#3288BD")
colors_agn     <- c("#2d6a4f", "orange3", "mediumpurple4")
colors_sf      <- c("Star-forming" = "#023e8a", "Quiescent" = "#d00000")
colors_morph   <- c("LTG" = "#023e8a", "ETG" = "#d00000")
colors_ma      <- c("#52b788", "#2d6a4f", "#081c15")
values         <- c("#d00000",  "#9d2730", "#3f37c9", "#052d5a")
cores          <- c("TN" = "#9E0142", "FN" = "#F46D43", "TP" = "#3288BD", "FP" = "#66C2A5")
cores          <- c("#9E0142", "#F46D43", "#3288BD", "#66C2A5")

# sizes: 
op_linewidth <- 1

# width and height figures:
width_figs  <- 12
height_figs <- 8

# variable names:
fSFG_label        <- expression(f[SFG])

label_logssfr       <- expression(log["10"]~(sSFR/yr^-1))
label_logsfr        <- expression(log["10"]~(SFR/M["☉"]~yr^-1))

label_logMstar      <- expression(log["10"]~(M["★"]/M["☉"]))
label_logMgroup     <- expression(log["10"]~(M["h"]/h^-1~M["☉"]))
label_logvelDisp_e  <- expression(log["10"]~(sigma/km~s^-1))
label_velDisp_e     <- expression(sigma~(km~s^-1))
label_logRproj_rvir <- expression(log["10"]~(R["proj"]/r["vir"]))
label_logSigma_SFR  <- expression(log["10"]~(Sigma["SFR"]))
label_TType         <- "T-Type"
label_conc          <- expression(R90[r]^{petro}/R50[r]^{petro})

label_rproj      <- expression(R["proj"]/r["vir"])
label_o3hb       <- expression("log"["10"] ~ "([OIII]"["\u03BB5007"] ~ "/H\u03B2)")
label_n2ha       <- expression("log"["10"] ~ "([NII]"["\u03BB6584"] ~ "/H\u03B1)")
label_eqw_ha     <- expression(""*EQW[H*alpha]*""~(ring(A)))
label_magPetro_r <- expression(Mag[Petro] * r)

# ggsave(path = paste0(figswd, "Mstar_sSFR_diagrams/Ma", Mhalo_min_group, "_zmax", z_max),
#        filename = paste0("Mstar_sSFR_", title, " - ", subtitle, ".pdf"),
#        device = cairo_pdf, width = width_figs, height = height_figs, units = "in", dpi = 600)
