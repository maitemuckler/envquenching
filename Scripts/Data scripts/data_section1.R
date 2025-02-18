## Pacotes ----
library(data.table)
library(ggplot2)
library(scales)
library(patchwork)
library(dplyr)
library(extrafont)
library(tidyr)

## Diretórios ----
wdcode <- "Scripts/"
wddata <- "/home/muckler/Work/Research/Astronomy/Data/"
wdfigs <- "Figures/"

## Outros códigos ----
source(file.path(wdcode, "Themes/my_theme.R"))
source(file.path(wdcode, "MyFunctions/sigma_gal_correction.R"))
source(file.path(wdcode, "MyFunctions/SF_Q_class.R"))
source(file.path(wdcode, "MyFunctions/distance_to_line_calc.R"))

## Configuração inicial ----
set.seed(123)
figs_width  <- 14
figs_height <- 8

## Definindo input e output files ----
input_file  <- "letter1_sample.csv"
output_file <- file.path(wddata, paste0("clean_", input_file))

## Lendo os dados ----
df <- fread(file.path(wddata, input_file))

## Tratamento de dados ----
# Substituindo valores inválidos por NA
replace_na_values <- function(data, columns, threshold) {
  data %>% mutate(across(all_of(columns), ~ if_else(.x <= threshold, NA_real_, .x)))
}

na_replacements <- list(
  list(columns = c("magPetro_u", "magPetro_g", "magPetro_r", "magPetro_i", "magPetro_z",
                   "magFiber_u", "magFiber_g", "magFiber_r", "magFiber_i", "magFiber_z",
                   "petroR90_r", "petroR50_r", "lgm_tot_p50", "lgm_tot_p16", "lgm_tot_p84", 
                   "sfr_tot_p50", "sfr_tot_p16", "sfr_tot_p84"), threshold = -9999),
  list(columns = c("logMstar", "logMstar_err", "logSFR_SED", "logSFR_SED_err"), threshold = -99)
)

for (replacement in na_replacements) {
  df <- replace_na_values(df, replacement$columns, replacement$threshold)
}

## Removendo NAs de colunas essenciais
df <- df %>% drop_na(absPetro_r, logMstar)

## Correção da dispersão de velocidades ----
df$velDisp_e    <- sigma_gal_correction(df$velDisp, df$Rhlr, df$Scale)
df$logvelDisp_e <- log10(df$velDisp_e)

# Removendo valores inválidos de sigma
df <- df %>% drop_na(velDisp_e) %>% filter(!is.infinite(logvelDisp_e))

# Corte em 30 km/s < velDisp_e < 320 km/s
ecdf_func <- ecdf(df$velDisp_e)
limits <- c(30, 320)
percentages <- round(ecdf_func(limits) * 100, 2)

# Gráfico de distribuição de velDisp_e
ggplot(subset(df, df$velDisp_e <= 500)) + 
  geom_histogram(aes(y = after_stat(density), x = velDisp_e), colour = 1, bins = 50, fill = "lightblue") + 
  geom_density(aes(velDisp_e), linewidth = 2, color = "grey20") +
  geom_vline(aes(xintercept = limits[1]), linetype = 2, linewidth = 2, color = "#F46D43") + 
  geom_vline(aes(xintercept = limits[2]), linetype = 2, linewidth = 2, color = "#F46D43") + 
  annotate("label", x = limits[1], y = 0, label = paste0("< ", limits[1], " km/s\n(", percentages[1], "% da amostra)"),
           color = "black", fill = "#F68765", size = 8, hjust = 0, vjust = -0.5) +
  annotate("label", x = limits[2], y = 0, label = paste0("< ", limits[2], " km/s\n(", percentages[2], "% da amostra)"),
           color = "black", fill = "#F68765", size = 8, hjust = 0, vjust = -0.5) +
  labs(x = expression(sigma(km~s^-1)), y = "Densidade") +
  theme_Publication()

df <- df %>% filter(between(velDisp_e, limits[1], limits[2]))

## Classificação SF/Q ----
df$logsSFR_GSWLC  <- df$logSFR_SED - df$logMstar
df$logsSFR_MPAJHU <- df$sfr_tot_p50 - df$lgm_tot_p50

mpajhu_params <- list(intercept = -7.85, slope = -0.3)
gswlc_params  <- list(intercept = -6.55, slope = -0.45)

df$SF_GSWLC  <- SF_Q_class(df$logsSFR_GSWLC,  df$logMstar, gswlc_params$slope, gswlc_params$intercept)
df$SF_MPAJHU <- SF_Q_class(df$logsSFR_MPAJHU, df$lgm_tot_p50, mpajhu_params$slope, mpajhu_params$intercept)

df$SF_GSWLC  <- factor(df$SF_GSWLC, levels = c("Star-forming", "Quiescent"))
df$SF_MPAJHU <- factor(df$SF_MPAJHU, levels = c("Star-forming", "Quiescent"))

## Criando variável dicotômica para T-Type ----
df$TType_label <- factor(ifelse(df$TType <= 0, "ETG", "LTG"))
df$morph_char  <- ifelse(df$TType_label == "ETG", "Early-type", "Late-type")

# Gráfico de distribuição de T-Type
ggplot(data = df, aes(x = TType)) + 
  geom_histogram(aes(y = after_stat(density), fill = TType_label), colour = "white", bins = 50) + 
  scale_fill_manual(values = colors_morph) +
  labs(x = "T-Type", y = "Densidade") +
  theme_Publication()

## Salvando os dados ----
write.csv(df, output_file, row.names = FALSE)
