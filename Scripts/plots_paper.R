library(ggplot2)
library(scales)
library(RColorBrewer)
library(wesanderson)
library(scico)
library(ggpubr)

source(paste0(wdmain, "Scripts/Themes/ggplot_theme_Publication-2.R"))

# Gradient color
pal <- wes_palette("Cavalcanti1", 10, type = "continuous")

# ssfr vs. stellar mass ----

gswlc_intercept  <- -6.55
gswlc_slope      <- -0.45

# ggplot(predictions, 
#        aes(x = logMstar, y = logsSFR_GSWLC)) + 
#   stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", contour = TRUE, bins = 30) + 
#   scale_fill_scico(palette = "devon", breaks = scales::pretty_breaks(n = 7), name = "Level") + 
#   scale_y_continuous(breaks = pretty_breaks(n = 10)) + 
#   scale_x_continuous(breaks = pretty_breaks(n = 10)) + 
#   #geom_abline(slope = -0.3, intercept = -7.85, color = 'red') + # Knobel+15
#   geom_abline(slope = gswlc_slope, intercept = gswlc_intercept, color = 'white', linetype = "dashed", linewidth = 2) + 
#   labs(x = mstar_label,
#        y = ssfr_label) +
#   theme_Publication() + 
#   theme(axis.text = element_text(size = 30),
#         axis.title = element_text(size = 36),
#         legend.position = "right",
#         legend.key.size = unit(1.5, "cm"),
#         legend.spacing = unit(3, "cm"),
#         legend.text = element_text(size = 24, margin = margin(l = 10)),
#         legend.title = element_text(size = 26, margin = margin(b = 15)),
#         legend.direction = "vertical",
#         legend.box = "horizontal")

ggplot(predictions, 
       aes(x = logMstar, y = logsSFR_GSWLC)) + 
  geom_point(aes(color = SF_char), size = 0.5) + 
  scale_color_manual(values = colors_sf) +
  scale_y_continuous(breaks = pretty_breaks(n = 10)) + 
  scale_x_continuous(breaks = pretty_breaks(n = 10)) + 
  geom_abline(slope = gswlc_slope, intercept = gswlc_intercept, color = 'black', linetype = "dashed", linewidth = 2) + 
  labs(x = mstar_label,
       y = ssfr_label) +
  theme_Publication() + 
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 36),
        legend.position = c(0.86, 0.88),
        legend.key.size = unit(1.5, "cm"),
        legend.spacing = unit(3, "cm"),
        legend.text = element_text(size = 24, margin = margin(l = 10)),
        legend.title = element_blank(),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.box.background = element_rect(colour = "black", linewidth = 2)) + 
  guides(color = guide_legend(override.aes = list(size = 4))) + 
  annotate(
    'text',
    x = 11.72,
    y = -10.72,
    label = 'Slope: -0.45\nIntercept: -6.55',
    fontface = 'italic', 
    size = 8,
    hjust = 0
  ) +
  annotate(
    'rect',
    xmin = 11.71,
    xmax = 12,
    ymin = -11.1,
    ymax = -10.36,
    alpha = 0, 
    col = 'black',
    linewidth = 1,
  ) +
  annotate(
    'curve',
    x = 11.7, # Play around with the coordinates until you're satisfied
    y = -10.7,
    xend = 11.6,
    yend = -11.5,
    linewidth = 2,
    curvature = 0.5,
    arrow = arrow(length = unit(0.5, 'cm'))
  )


# --------------------------------------------------------------------------------------------------------------------------------------------

df_means <- aggregate(e ~ logMstar_cut + logvelDisp_cut + misclass1, data = dados, mean)

g5 <- ggplot(df_means, aes(x = logMstar_cut, y = logvelDisp_cut, fill = e)) +
  geom_tile() + 
  scale_fill_viridis(option = "inferno", direction = -1) + 
  facet_wrap( . ~ misclass1) + 
  scale_x_discrete(breaks = logMstar_midpoints[seq(1, length(logMstar_midpoints), by = 4)], 
                   labels = round(logMstar_midpoints[seq(1, length(logMstar_midpoints), by = 4)], digits = 2)) + 
  scale_y_discrete(breaks = logvelDisp_midpoints[seq(1, length(logvelDisp_midpoints), by = 4)],
                   labels = round(logvelDisp_midpoints[seq(1, length(logvelDisp_midpoints), by = 4)], digits = 2)) + 
  theme_bw()

g5

ggsave(plot = g5, path = wdfigs, file = 'g5.pdf', 
       width = 12, height = 8, units = "cm", device = cairo_pdf)

df_means <- aggregate(B_T_r ~ logMstar_cut + logvelDisp_cut + misclass1, data = dados, mean)

g5 <- ggplot(df_means, aes(x = logMstar_cut, y = logvelDisp_cut, fill = B_T_r)) +
  geom_tile() + 
  scale_fill_viridis(option = "inferno", direction = -1) + 
  facet_wrap( . ~ misclass1) + 
  scale_x_discrete(breaks = logMstar_midpoints[seq(1, length(logMstar_midpoints), by = 4)], 
                   labels = round(logMstar_midpoints[seq(1, length(logMstar_midpoints), by = 4)], digits = 2)) + 
  scale_y_discrete(breaks = logvelDisp_midpoints[seq(1, length(logvelDisp_midpoints), by = 4)],
                   labels = round(logvelDisp_midpoints[seq(1, length(logvelDisp_midpoints), by = 4)], digits = 2)) + 
  theme_bw()

g5

df_means <- aggregate(Nth_gal_ra_Neq5 ~ logMstar_cut + logvelDisp_cut + misclass1, data = dados, mean)

g5 <- ggplot(df_means, aes(x = logMstar_cut, y = logvelDisp_cut, fill = Nth_gal_ra_Neq5)) +
  geom_tile() + 
  scale_fill_viridis(option = "inferno", direction = -1) + 
  facet_wrap( . ~ misclass1) + 
  scale_x_discrete(breaks = logMstar_midpoints[seq(1, length(logMstar_midpoints), by = 4)], 
                   labels = round(logMstar_midpoints[seq(1, length(logMstar_midpoints), by = 4)], digits = 2)) + 
  scale_y_discrete(breaks = logvelDisp_midpoints[seq(1, length(logvelDisp_midpoints), by = 4)],
                   labels = round(logvelDisp_midpoints[seq(1, length(logvelDisp_midpoints), by = 4)], digits = 2)) + 
  theme_bw()

g5

df_means <- aggregate(logRproj_rvir ~ logMstar_cut + logvelDisp_cut + misclass1, data = dados, mean)

g5 <- ggplot(df_means, aes(x = logMstar_cut, y = logvelDisp_cut, fill = logRproj_rvir)) +
  geom_tile() + 
  scale_fill_viridis(option = "inferno", direction = -1) + 
  facet_wrap( . ~ misclass1) + 
  scale_x_discrete(breaks = logMstar_midpoints[seq(1, length(logMstar_midpoints), by = 4)], 
                   labels = round(logMstar_midpoints[seq(1, length(logMstar_midpoints), by = 4)], digits = 2)) + 
  scale_y_discrete(breaks = logvelDisp_midpoints[seq(1, length(logvelDisp_midpoints), by = 4)],
                   labels = round(logvelDisp_midpoints[seq(1, length(logvelDisp_midpoints), by = 4)], digits = 2)) + 
  theme_bw()

g5

df_means <- aggregate(P_disk ~ logMstar_cut + logvelDisp_cut + misclass1, data = dados, mean)

g6 <- ggplot(df_means, aes(x = logMstar_cut, y = logvelDisp_cut, fill = P_disk)) +
  geom_tile() + 
  scale_fill_viridis(option = "inferno", direction = -1) + 
  facet_wrap( . ~ misclass1) + 
  scale_x_discrete(breaks = logMstar_midpoints[seq(1, length(logMstar_midpoints), by = 4)], 
                   labels = round(logMstar_midpoints[seq(1, length(logMstar_midpoints), by = 4)], digits = 2)) + 
  scale_y_discrete(breaks = logvelDisp_midpoints[seq(1, length(logvelDisp_midpoints), by = 4)],
                   labels = round(logvelDisp_midpoints[seq(1, length(logvelDisp_midpoints), by = 4)], digits = 2)) + 
  theme_bw()

g6

ggsave(plot = g6, path = wdfigs, file = 'g6.pdf', 
       width = 12, height = 8, units = "cm", device = cairo_pdf)

df_means <- aggregate(distLine_GSWLC ~ logMstar_cut + logvelDisp_cut + misclass1, data = dados, mean)

g6 <- ggplot(df_means, aes(x = logMstar_cut, y = logvelDisp_cut, fill = abs(distLine_GSWLC))) +
  geom_tile() + 
  scale_fill_viridis(option = "inferno", direction = -1) + 
  facet_wrap( . ~ misclass1) + 
  scale_x_discrete(breaks = logMstar_midpoints[seq(1, length(logMstar_midpoints), by = 4)], 
                   labels = round(logMstar_midpoints[seq(1, length(logMstar_midpoints), by = 4)], digits = 2)) + 
  scale_y_discrete(breaks = logvelDisp_midpoints[seq(1, length(logvelDisp_midpoints), by = 4)],
                   labels = round(logvelDisp_midpoints[seq(1, length(logvelDisp_midpoints), by = 4)], digits = 2)) + 
  theme_bw()

g6

ggsave(plot = g6, path = wdfigs, file = 'g6.pdf', 
       width = 12, height = 8, units = "cm", device = cairo_pdf)


df_means1 <- aggregate(logSFR_SED ~ logMstar_cut + logvelDisp_cut + misclass1, data = dados, mean)
df_means2 <- aggregate(sfr_tot_p50 ~ logMstar_cut + logvelDisp_cut + misclass1, data = dados, mean)

colnames(df_means1)[which(colnames(df_means1) == "logSFR_SED")]  <- 'logSFR'
colnames(df_means2)[which(colnames(df_means2) == "sfr_tot_p50")] <- 'logSFR'

df_means1$catalog <- "GSWLC"
df_means2$catalog <- "MPA-JHU"
df_means <- rbind(df_means1, df_means2)

g7 <- ggplot(df_means, aes(x = logMstar_cut, y = logvelDisp_cut, fill = logSFR)) +
  geom_tile() + 
  scale_fill_viridis(option = "inferno", direction = -1, breaks = pretty_breaks(n = 10)) + 
  facet_grid(catalog ~ misclass1) + 
  scale_x_discrete(breaks = logMstar_midpoints[seq(1, length(logMstar_midpoints), by = 8)], 
                   labels = round(logMstar_midpoints[seq(1, length(logMstar_midpoints), by = 8)], digits = 2)) + 
  scale_y_discrete(breaks = logvelDisp_midpoints[seq(1, length(logvelDisp_midpoints), by = 8)],
                   labels = round(logvelDisp_midpoints[seq(1, length(logvelDisp_midpoints), by = 8)], digits = 2)) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm"),            # Tamanho dos quadrados/círculos na legenda
        legend.spacing = unit(0.5, "cm"))

g7

ggsave(plot = g7, path = wdfigs, file = 'g7.pdf', 
       width = 20, height = 11, units = "cm", device = cairo_pdf)
