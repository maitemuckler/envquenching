## Lendo pacotes ----
library(data.table)
library(ggplot2)
library(scales)
library(patchwork)
library(dplyr)
library(extrafont)

# Carregue as fontes
#font_import()
#loadfonts(device = "postscript")

figs_width  <- 14
figs_height <- 8

## Diretórios ----
wdcode <- "Scripts/"
wddata <- "/home/muckler/Work/Research/Astronomy/Data/"
wdfigs <- "Figures/"

## Outros códigos ----
source(paste0(wdcode, "Themes/my_theme.R"))
source(paste0(wdcode, "MyFunctions/sigma_gal_correction.R"))
source(paste0(wdcode, "MyFunctions/SF_Q_class.R"))
source(paste0(wdcode, "MyFunctions/distance_to_line_calc.R"))

## Seed ----
set.seed(123)

## Definindo input e output files ----
input_file  <- "letter1_sample.csv"
output_file <- paste0("clean_", input_file)

## Lendo os dados ----
df <- fread(paste0(wddata, input_file)) # 287,248 galaxies

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
if(any(is.infinite(df$logvelDisp_e))){df <- df[-which(is.infinite(df$logvelDisp_e)),]} # 270.404

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

ggsave(path = paste0(wdfigs, "pdf/"),
       filename = "logMstar_logsSFR.pdf",
       device = cairo_pdf, 
       width = figs_width, height = figs_height, 
       units = "in", dpi = 600)

ggsave(path = paste0(wdfigs, "eps/"),
       filename = "logMstar_logsSFR.eps",
       device = cairo_ps,
       width = figs_width, height = figs_height,
       units = "in", dpi = 600)

ggsave(path = paste0(wdfigs, "png/"),
       filename = paste0("cut_velDisp.png"),
       width = figs_width, height = figs_height, 
       units = "in", dpi = 600)

df <- df %>%
  filter(velDisp_e > 30 & velDisp_e < 320) # 255.727 (- 14.677 galaxies)

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
  labs(x = label_logMstar,
       y = label_logssfr) +
  theme_Publication() + 
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 36),
        legend.position = "inside", 
        legend.position.inside =  c(0.15, 0.15),
        legend.key.size = unit(1.5, "cm"),
        legend.spacing = unit(3, "cm"),
        legend.text = element_text(size = 24, margin = margin(l = 10)),
        legend.title = element_blank(),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.box.background = element_rect(colour = "black", linewidth = 2)) + 
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) 

ggsave(path = paste0(wdfigs, "pdf/"),
       filename = "logMstar_logsSFR.pdf",
       device = cairo_pdf, 
       width = figs_width, height = figs_height,
       units = "in", dpi = 600)

ggsave(path = paste0(wdfigs, "eps/"),
       filename = "logMstar_logsSFR.eps",
       device = cairo_ps,
       width = figs_width, height = figs_height,
       units = "in", dpi = 600)

ggsave(path = paste0(wdfigs, "png/"),
       filename = "logMstar_logsSFR.png",
       width = figs_width, height = figs_height,
       units = "in", dpi = 600)

# Verificando valores NA 
names(df)[sapply(df, function(x) any(is.na(x)))]

# Verificando valores infinito
names(df)[sapply(df, function(x) any(is.infinite(x)))]

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
        legend.position = "inside", 
        legend.position.inside = c(0.8, 0.85),
        legend.key.size = unit(1.5, "cm"),
        legend.spacing = unit(3, "cm"),
        legend.text = element_text(size = 24, margin = margin(l = 10)),
        legend.title = element_blank(),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.box.background = element_rect(colour = "black", linewidth = 2)) + 
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) 

ggsave(path = paste0(wdfigs, "pdf/"),
       filename = "TType.pdf",
       device = cairo_pdf, 
       width = figs_width, height = figs_height,
       units = "in", dpi = 600)

ggsave(path = paste0(wdfigs, "eps/"),
       filename = "TType.eps",
       device = cairo_ps,
       width = figs_width, height = figs_height,
       units = "in", dpi = 600)

ggsave(path = paste0(wdfigs, "png/"),
       filename = "TTYpe.png",
       width = figs_width, height = figs_height,
       units = "in", dpi = 600)

## Salvar tabela ----
write.csv(df, paste0(wddata, output_file), row.names = F) # 255,727 galaxies


