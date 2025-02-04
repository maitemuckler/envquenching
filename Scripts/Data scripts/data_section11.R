## Classificação AGN ----
df$logO3Hb       <- log10(df$oiii_5007_flux/df$h_beta_flux)
df$logN2Ha       <- log10(df$nii_6584_flux/df$h_alpha_flux)

df$h_alpha_eqw  <- -df$h_alpha_eqw
df$nii_6548_eqw <- -df$nii_6548_eqw

aa <- 0.47; bb <- 1.19 # BPT (Kewley+01)

df$AGN <- df$logO3Hb <= (0.61/(df$logN2Ha - aa)) + bb &
  df$logN2Ha < 0.47

df$AGN[which(df$AGN == TRUE)]  <- "Non-AGN"
df$AGN[which(df$AGN == FALSE)] <- "AGN"

gal_nonAGN <- df[which(df$AGN == "Non-AGN"),]
gal_AGN    <- df[which(df$AGN == "AGN"),]
gal_NA     <- df[which(is.na(df$AGN)),]

# WHAN (p/ AGN no BPT)
gal_AGN$AGN[gal_AGN$h_alpha_eqw >= 3 & gal_AGN$logN2Ha >= -0.4] <- "AGN"
gal_AGN$AGN[gal_AGN$h_alpha_eqw < 3] <- "Non-AGN"
gal_AGN$AGN[gal_AGN$logN2Ha < -0.4]  <- "Non-AGN"

# WHAN(p/ NAs no BPT)
gal_NA$AGN[gal_NA$h_alpha_eqw >= 3 & gal_NA$logN2Ha >= -0.4] <- "AGN"
gal_NA$AGN[gal_NA$h_alpha_eqw < 3] <- "Non-AGN"
gal_NA$AGN[gal_NA$logN2Ha < -0.4]  <- "Non-AGN"

# Une as 3 tabelas novamente
df <- rbind(gal_nonAGN, gal_AGN, gal_NA)
rm(gal_AGN, gal_NA, gal_nonAGN)

# Removendo AGN NA verdadeiro

agn_na <- df[which(is.na(df$AGN)),]
nrow(agn_na)/nrow(df)*100

agn_na <- agn_na %>%
  select(AGN, logO3Hb, logN2Ha, oiii_5007_flux, h_beta_flux, nii_6584_flux, h_alpha_flux, h_alpha_eqw, nii_6548_eqw)

names(df)[sapply(df, function(x) any(is.infinite(x)))]
names(df)[sapply(df, function(x) any(is.na(x)))]

if(any(is.na(df$AGN))){df <- df[-which(is.na(df$AGN)),]}

df$AGN <- factor(df$AGN)
table(df$AGN)
(table(df$AGN)/nrow(df)) * 100
0.25*nrow(df) #63974

bpt <- ggplot(df[sample(1:nrow(df), size = 0.25*nrow(df)),], 
              aes(x = logN2Ha, y = logO3Hb, color = AGN)) + 
  geom_point(aes(alpha = AGN, size = AGN), na.rm = TRUE) + 
  
  scale_x_continuous(limits = c(-1.4, 0.6), breaks = pretty_breaks(n = 10)) +
  scale_y_continuous(limits = c(-2, 2), breaks = pretty_breaks(n = 10)) +
  
  scale_color_manual(values = c("#F46D43", "grey40"), name = "") +
  scale_alpha_manual(values = c(1, 0.1), name = "") + 
  scale_size_manual(values = c(1, 0.1), name = "") + 
  labs(y = o3hb_label) + 
  theme_Publication() +
  theme(axis.text.y = element_text(size = 30),
        axis.title.y = element_text(size = 32, margin = margin(r = 15)),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #legend.position.inside = c(0.16, 0.91),
        legend.position = "top",
        legend.key.size = unit(1.5, "cm"),
        legend.spacing = unit(3, "cm"),
        legend.text = element_text(size = 24, margin = margin(l = 10)),
        legend.title = element_text(size = 26, margin = margin(b = 15)),
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.box.background = element_rect(colour = "black", linewidth = 2),
        plot.margin = margin(t = 10, r = 10, b = 0, l = 10)) + 
  guides(color = guide_legend(override.aes = list(size = 6, alpha = 1))) + 
  
  annotate(
    'text',
    x = 0.5,
    y = 1.8,
    label = 'BPT diagram',
    fontface = 'italic', 
    size = 8
  )

# Gráfico WHAN
whan <- ggplot(df[sample(1:nrow(df), size = 0.25*nrow(df)),], 
               aes(x = logN2Ha, y = h_alpha_eqw, color = AGN)) + 
  geom_point(aes(alpha = AGN, size = AGN), na.rm = TRUE) + 
  
  scale_x_continuous(limits = c(-1.4, 0.6), breaks = pretty_breaks(n = 10)) + 
  scale_y_continuous(limits = c(0.1, 100), breaks = c(0.5, 1, 10, 5, 50, 100)) +
  
  scale_color_manual(values = c("#F46D43", "grey40"), name = "") +
  scale_alpha_manual(values = c(1, 0.1), name = "") + 
  scale_size_manual(values = c(1, 0.1), name = "") + 
  coord_trans(y = "log10") +
  labs(y = eqw_ha_label, x = n2ha_label) + 
  theme_Publication() +
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 32),
        axis.title.y = element_text(margin = margin(r = 5)), # Ajusta a margem do título do eixo y
        legend.position = "none",
        plot.margin = margin(t = 0, r = 10, b = 10, l = 10)) +
  
  annotate(
    'text',
    x = 0.5,
    y = 80,
    label = 'WHAN diagram',
    fontface = 'italic', 
    size = 8
  )

agn <- bpt / whan
agn

ggsave(path = wdfigs,
       filename = paste0("AGN.pdf"),
       device = cairo_pdf, width = 14, height = 16, 
       units = "in", dpi = 600)

ggsave(path = paste0(wdfigs, "png/"),
       filename = paste0("AGN.png"),
       width = 16, height = 10, 
       units = "in", dpi = 600)

## Renomear e criar variáveis ----

# ----------- Renomear:
colnames(df)[which(colnames(df) == "__B_T_r")] <- "B_T_r"
colnames(df)[which(colnames(df) == "e__B_T_r")] <- "e_B_T_r"

# ----------- Criar: 
## Criar distância até a reta 
df$distLine_GSWLC  <- distance_to_line(x = df$logMstar, y = df$logsSFR_GSWLC, slope = gswlc_slope, intercept = gswlc_intercept)
df$distLine_MPAJHU <- distance_to_line(x = df$lgm_tot_p50, y = df$logsSFR_MPAJHU, slope = mpajhu_slope, intercept = mpajhu_intercept)

# Example:
# mstar = 10.5
# ssfr = -10
# numerador   = - ((0.45*mstar) + ssfr + 6.55)
# numerador2  = ((-0.45*mstar) - ssfr - 6.55)
# denominador = sqrt((-0.45)^2 + 1)
# denominador2 = sqrt((0.45)^2 + 1)
# numerador/denominador
# numerador2/denominador
# numerador2/denominador2

## Criar SFR Surface Density 
# https://arxiv.org/pdf/2310.11493
logR_iso <- (0.188*df$logMstar) + (0.333*log10(df$Rhlr)) - 1.037 # (eq. 6)
R_iso    <- 10^logR_iso
A       <- 1.26*pi*R_iso^2

df$Sigma_SFR    <- (10^df$logSFR_SED)/A
df$logSigma_SFR <- log10(df$Sigma_SFR)

ggplot(data = df[sample(1:nrow(df), size = 0.25*nrow(df)),], 
       aes(x = logMstar, y = logSigma_SFR)) + 
  geom_point(aes(color = SF_GSWLC), size = 0.5, alpha = 0.1) + 
  scale_color_manual(values = colors_sf) +
  scale_y_continuous(breaks = pretty_breaks(n = 10)) + 
  scale_x_continuous(breaks = pretty_breaks(n = 10)) + 
  geom_abline(slope = gswlc_slope, intercept = gswlc_intercept, color = 'black', linetype = "dashed", linewidth = 2) + 
  labs(x = label_logMstar,
       y = label_logSigma_SFR) +
  theme_Publication() + 
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 36),
        legend.position.inside = c(0.15, 0.15),
        legend.key.size = unit(1.5, "cm"),
        legend.spacing = unit(3, "cm"),
        legend.text = element_text(size = 24, margin = margin(l = 10)),
        legend.title = element_blank(),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.box.background = element_rect(colour = "black", linewidth = 2)) + 
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) 

ggsave(path = wdfigs,
       filename = paste0("Sigma_SFR.pdf"),
       device = cairo_pdf, width = 16, height = 10, 
       units = "in", dpi = 600)

ggsave(path = paste0(wdfigs, "png/"),
       filename = paste0("Sigma_SFR.png"),
       width = 16, height = 10, 
       units = "in", dpi = 600)

## Criar concentração 
df$conc <- df$petroR90_r/df$petroR50_r
