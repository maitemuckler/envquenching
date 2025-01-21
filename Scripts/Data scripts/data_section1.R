## Lendo pacotes ----
library(data.table)
library(ggplot2)
library(scales)
library(patchwork)

## Diretórios ----
wdcode <- "~/Work/Research/Astronomy/Projects/environmental-quenching/Scripts/"
wddata <- "~/Work/Research/Astronomy/Data/"
wdfigs <- "~/Work/Research/Astronomy/Projects/environmental-quenching/Figures/"

## Outros códigos ----
source(paste0(wdcode, "Themes/my_theme.R"))
source(paste0(wdcode, "MyFunctions/sigma_gal_correction.R"))
source(paste0(wdcode, "MyFunctions/SF_Q_class.R"))
source(paste0(wdcode, "MyFunctions/distance_to_line_calc.R"))

## Seed ----
set.seed(123)

## Definindo input e output files ----
input_file  <- "SDSS_DR18_Legacy_MGS_QSO_localDensity+GSWLC+Simard11+DS18.csv"
output_file <- "clean_SDSS_DR18_Legacy_MGS_QSO_localDensity+GSWLC+Simard11+DS18.csv"

## Lendo os dados ----
df <- fread(paste0(wddata, "EnvQuenching/", input_file)) # 287,346 galaxies

## Substituindo -9999 por NA (sem remover linhas) ----

colunas_para_substituir <- c("magPetro_u", "magPetro_g", "magPetro_r", "magPetro_i", "magPetro_z",
                             "magFiber_u", "magFiber_g", "magFiber_r", "magFiber_i", "magFiber_z",
                             "petroR90_r", "petroR50_r", 
                             "lgm_tot_p50", "lgm_tot_p16", "lgm_tot_p84", 
                             "sfr_tot_p50", "sfr_tot_p16", "sfr_tot_p84")

df <- df %>%
  mutate(across(all_of(colunas_para_substituir), ~ if_else(.x <= -9999, NA_real_, .x)))

# Susbtituir valores -99 por NA:
colunas_para_substituir <- c("logMstar", "logMstar_err", "logSFR_SED", "logSFR_SED_err")
df <- df %>%
  mutate(across(all_of(colunas_para_substituir), ~ if_else(.x <= -99, NA_real_, .x)))

rm(colunas_para_substituir)
summary(df)

## Limpeza e correção de dispersão de velocidades ----

names(df)[sapply(df, function(x) any(is.infinite(x)))] # Nenhum infinito
names(df)[sapply(df, function(x) any(is.na(x)))]

# Removendo NA:
if(any(is.na(df$absPetro_r))){df <- df[-which(is.na(df$absPetro_r)),]} # 287.326
if(any(is.na(df$logMstar))){df <- df[-which(is.na(df$logMstar)),]} # 285.709

# Corrigindo sigma:
df$velDisp_e    <- sigma_gal_correction(df$velDisp, df$Rhlr, df$Scale)
df$logvelDisp_e <- log10(df$velDisp_e)

names(df)[sapply(df, function(x) any(is.infinite(x)))]
names(df)[sapply(df, function(x) any(is.na(x)))]

# Removendo Na e infinitos de sigma:
if(any(is.na(df$velDisp_e))){df <- df[-which(is.na(df$velDisp_e)),]} # 285.692
if(any(is.infinite(df$logvelDisp_e))){df <- df[-which(is.infinite(df$logvelDisp_e)),]} # 270.507

# Corte em 30 km/s < velDisp_e < 320 km/s:
ecdf_func <- ecdf(df$velDisp_e) # Calcular a função de distribuição acumulada empírica

sigma_30      <- 30 # Valor específico para o qual queremos saber a porcentagem da amostra
sigma_30_perc <- ecdf_func(sigma_30) * 100 # Calcular a porcentagem da amostra que é menor ou igual ao valor específico
sigma_30_perc <- round(sigma_30_perc, digits = 2)

sigma_320      <- 320 
sigma_320_perc <- ecdf_func(sigma_320) * 100 
sigma_320_perc <- round(sigma_320_perc, digits = 2)

ggplot(subset(df, df$velDisp_e <= 500)) + 
  geom_histogram(aes(y = after_stat(density), x = velDisp_e), colour = 1, bins = 50, fill = "lightblue") + 
  geom_density(aes(velDisp_e), linewidth = 2, color = "grey20") +
  geom_vline(aes(xintercept = 30), linetype = 2, linewidth = 2, color = "#F46D43") + 
  geom_vline(aes(xintercept = 320), linetype = 2, linewidth = 2, color = "#F46D43") + 
  scale_x_continuous(breaks = pretty_breaks(n = 10)) + 
  scale_y_continuous(breaks = pretty_breaks(n = 10)) + 
  
  annotate("label", x = sigma_30, y = 0, 
           label = paste0("< 30 km/s\n(",sigma_30_perc,"% of the sample)"),
           color = "black", fill = "#F68765", size = 8, hjust = 0, vjust = -0.5, label.padding = unit(0.25, "lines")) +
  
  annotate("label", x = sigma_320, y = 0, 
           label = paste0("< 320 km/s\n(",sigma_320_perc,"% of the sample)"),
           color = "black", fill = "#F68765", size = 8, hjust = 0, vjust = -0.5, label.padding = unit(0.25, "lines")) +
  
  labs(x = expression(sigma(km~s^-1)),
       y = "Density") +
  theme_Publication() + 
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 36, face = "plain"),
        legend.position.inside = c(0.15, 0.15),
        legend.key.size = unit(1.5, "cm"),
        legend.spacing = unit(3, "cm"),
        legend.text = element_text(size = 24, margin = margin(l = 10)),
        legend.title = element_blank(),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.box.background = element_rect(colour = "black", linewidth = 2))

ggsave(path = wdfigs,
       filename = paste0("cut_velDisp.pdf"),
       device = cairo_pdf, width = 16, height = 10, 
       units = "in", dpi = 600)

ggsave(path = paste0(wdfigs, "png/"),
       filename = paste0("cut_velDisp.png"),
       width = 16, height = 10, 
       units = "in", dpi = 600)

df <- df %>%
  filter(velDisp_e > 30 & velDisp_e < 320) # 255.825 (- 14.682 galaxies)

## Classificação SF/Q ----

df$logsSFR_GSWLC  <- df$logSFR_SED - df$logMstar
df$logsSFR_MPAJHU <- df$sfr_tot_p50 - df$lgm_tot_p50

mpajhu_intercept <- -7.85
mpajhu_slope     <- -0.3
gswlc_intercept  <- -6.55
gswlc_slope      <- -0.45

df$SF_GSWLC  <- SF_Q_class(df$logsSFR_GSWLC,  df$logMstar, gswlc_slope, gswlc_intercept)
df$SF_MPAJHU <- SF_Q_class(df$logsSFR_MPAJHU, df$lgm_tot_p50, mpajhu_slope, mpajhu_intercept)

df$SF_GSWLC  <- factor(df$SF_GSWLC, levels = c("Star-forming", "Quiescent"))
df$SF_MPAJHU <- factor(df$SF_MPAJHU, levels = c("Star-forming", "Quiescent"))

table(df$SF_GSWLC)
(table(df$SF_GSWLC)/nrow(df)) * 100
0.25*nrow(df)

ggplot(df[sample(1:nrow(df), size = 0.25*nrow(df)),], 
       aes(x = logMstar, y = logsSFR_GSWLC)) + 
  geom_point(aes(color = SF_GSWLC), size = 0.5, alpha = 0.1) + 
  scale_color_manual(values = colors_sf) +
  scale_y_continuous(breaks = pretty_breaks(n = 10)) + 
  scale_x_continuous(breaks = pretty_breaks(n = 10)) + 
  geom_abline(slope = gswlc_slope, intercept = gswlc_intercept, color = 'black', linetype = "dashed", linewidth = 2) + 
  labs(x = mstar_label,
       y = ssfr_label) +
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
       filename = paste0("logMstar_logsSFR.pdf"),
       device = cairo_pdf, width = 16, height = 10, 
       units = "in", dpi = 600)

ggsave(path = paste0(wdfigs, "png/"),
       filename = paste0("logMstar_logsSFR.png"),
       width = 16, height = 10, 
       units = "in", dpi = 600)

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

# Verificando valores NA 
names(df)[sapply(df, function(x) any(is.na(x)))]

# Verificando valores infinito
names(df)[sapply(df, function(x) any(is.infinite(x)))]

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

## Criar T-Type dicotômico

# T-Type < 0 correspond to ETGs,
# T-Type = 0 are S0, 
# T-Type > 0 are spiral galaxies (from Sa to Sm), 
# T-Type=10 are irregular galaxies.

df$TType_label <- NA
df$TType_label[which(df$TType <= 0)] <- "ETG"
df$TType_label[which(df$TType > 0)]  <- "LTG"
df$TType_label <- as.factor(df$TType_label)

ggplot(data = df, aes(x = TType)) + 
  geom_histogram(aes(y = after_stat(density), fill = TType_label), colour = "white", bins = 50) + 
  scale_fill_manual(values = colors_morph, labels = c("ETG" = "Early-type galaxies (ETG)",
                                                      "LTG" = "Late-type galaxies (LTG)")) +
  scale_y_continuous(breaks = pretty_breaks(n = 10)) + 
  scale_x_continuous(breaks = pretty_breaks(n = 10)) + 
  labs(x = "T-Type",
       y = "Density") +
  theme_Publication() + 
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 36),
        legend.position.inside = c(0.8, 0.85),
        legend.key.size = unit(1.5, "cm"),
        legend.spacing = unit(3, "cm"),
        legend.text = element_text(size = 24, margin = margin(l = 10)),
        legend.title = element_blank(),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.box.background = element_rect(colour = "black", linewidth = 2)) + 
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) 

ggsave(path = wdfigs,
       filename = paste0("TType.pdf"),
       device = cairo_pdf, width = 16, height = 10, 
       units = "in", dpi = 600)

ggsave(path = paste0(wdfigs, "png/"),
       filename = paste0("TTYpe.png"),
       width = 16, height = 10, 
       units = "in", dpi = 600)

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
  labs(x = mstar_label,
       y = Sigma_SFR_label) +
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

## Salvar tabela ----
write.csv(df, paste0(wddata, "EnvQuenching/", output_file), row.names = F) # 287,346 galaxies