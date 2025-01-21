setwd("~/Work/Research/Astronomy/Projects/environmental-quenching/")

library(data.table)
library(binom)
library(ggplot2)
library(InformationValue)
library(caret)
library(dplyr)
library(scales)
library(viridis)
library(ggthemes) # theme_clean

source("Scripts/Themes/my_theme.R")

wddata <- "~/Work/Research/Astronomy/Data/environmental-quenching-data/"

zmax          <- 0.03
catalog       <- "GSWLC"
TType_lim     <- 0
critval       <- 1.96 
Ma_vec <- c(12.3, 13, 14)

SF_name      <- ifelse(catalog == "GSWLC", "SF_GSWLC", "SF_MPAJHU")
logsSFR_name <- ifelse(catalog == "GSWLC", "logsSFR_GSWLC", "logsSFR_MPAJHU")

df <- data.frame()
columns_names <- c("groupID", "igal", "ra", "dec", "z", SF_name, logsSFR_name, "TType", 
                   "type", "logvelDisp_e", 
                   "logRproj_rvir", "logMgroup")

for (i in 1:3) {
  Ma         <- Ma_vec[i]
  input_data <- paste0("/inputdata_zmax", zmax, "_Rlim2.5_Ma", Ma, "_flag_good==1_MANGLE_logMstar_min10.5.csv")
  aux        <- fread(paste0(wddata, "inputModel/", catalog, input_data), select = columns_names)
  aux$Ma     <- Ma
  df         <- rbind(df, aux)
}

rm(aux)

colnames(df)[which(colnames(df) == SF_name)] <- "SF_char"
df$SF_char  <- as.factor(df$SF_char)
df$SF       <- ifelse(df$SF_char == "Star-forming", 1, 0)

colnames(df)[which(colnames(df) == logsSFR_name)] <- "logsSFR"

df$LT <- ifelse(df$TType >= TType_lim, 1, 0)
df$LT <- as.factor(df$LT)

df$velDisp_e <- 10^df$logvelDisp_e

data.s <- subset(df, df$type == "Satellite")
data.c <- subset(df, df$type == "Central")

# Defina aqui os limites dos bins de sigma
sigma_breaks <- c(min(data.s$velDisp_e), 130, 150, max(data.s$velDisp_e))
sigma_breaks <- c(min(data.s$velDisp_e), 60, 100, max(data.s$velDisp_e))  

data.s_cut      <- subset(data.s, data.s$velDisp_e < 130 & data.s$velDisp_e > 70)
logRproj_breaks <- quantile(data.s_cut$logRproj_rvir, probs = seq(0,1, by = 0.1))

data_binned.sSFR <- data.s_cut %>%
  mutate(
    logRproj_rvir_bin = cut(logRproj_rvir, breaks = logRproj_breaks, include.lowest = TRUE),
    logMgroup_bin = cut(logMgroup, breaks = c(12.3, 13.5, 14.0, 15.5), include.lowest = TRUE)) %>%
  group_by(logRproj_rvir_bin, logMgroup_bin) %>%
  summarise(
    mean_logRproj_rvir = mean(logRproj_rvir, na.rm = TRUE), 
    mean_logsSFR = mean(logsSFR, na.rm = TRUE),
    sd_logsSFR = sd(logsSFR, na.rm = TRUE),  # Desvio padrão
    n = n()  # Número de observações em cada bin
  ) %>%
  mutate(
    error_logsSFR = sd_logsSFR / sqrt(n)  # Erro padrão da média
  ) %>%
  ungroup()

data_binned.TType <- data.s_cut %>%
  mutate(
    logRproj_rvir_bin = cut(logRproj_rvir, breaks = logRproj_breaks, include.lowest = TRUE),
    logMgroup_bin = cut(logMgroup, breaks = c(12.3, 13.5, 14.0, 15.5), include.lowest = TRUE)) %>%
  group_by(logRproj_rvir_bin, logMgroup_bin) %>%
  summarise(
    mean_logRproj_rvir = mean(logRproj_rvir, na.rm = TRUE), 
    mean_TType = mean(TType, na.rm = TRUE),
    sd_TType = sd(TType, na.rm = TRUE),  # Desvio padrão
    n = n()  # Número de observações em cada bin
  ) %>%
  mutate(
    error_TType = sd_TType / sqrt(n)  # Erro padrão da média
  ) %>%
  ungroup()

data_binned.velDisp_e <- data.s_cut %>%
  mutate(
    logRproj_rvir_bin = cut(logRproj_rvir, breaks = logRproj_breaks, include.lowest = TRUE),
    logMgroup_bin = cut(logMgroup, breaks = c(12.3, 13.5, 14.0, 15.5), include.lowest = TRUE)) %>%
  group_by(logRproj_rvir_bin, logMgroup_bin) %>%
  summarise(
    mean_logRproj_rvir = mean(logRproj_rvir, na.rm = TRUE), 
    mean_velDisp_e = mean(velDisp_e, na.rm = TRUE),
    sd_velDisp_e = sd(velDisp_e, na.rm = TRUE),  # Desvio padrão
    n = n()  # Número de observações em cada bin
  ) %>%
  mutate(
    error_velDisp_e = sd_velDisp_e / sqrt(n)  # Erro padrão da média
  ) %>%
  ungroup()

# Criação do gráfico com barras de erro
ggplot(data_binned.sSFR, aes(x = mean_logRproj_rvir, y = mean_logsSFR, color = logMgroup_bin)) +
  geom_line(aes(group = logMgroup_bin), size = 0.8) +  # Linhas conectando os pontos
  geom_point(shape = 15, size = 3) +  # Pontos quadrados
  geom_errorbar(aes(ymin = mean_logsSFR - error_logsSFR, ymax = mean_logsSFR + error_logsSFR),
                width = 0.05, size = 0.5) +  # Barras de erro
  scale_color_manual(values = c("red", "orange", "blue"),  
                     labels = c("12.3 < log Mh < 13.5",
                                "13.5 < log Mh < 14.0",
                                "14.0 < log Mh < 15.5")) +
  scale_y_continuous(breaks = pretty_breaks(n = 6)) + 
  scale_x_continuous(breaks = pretty_breaks(n = 6)) + 
  labs(x = expression(log(R[proj] / r[vir])),
       y = expression("Mean log sSFR"),
       color = NULL) +  # Remove título da legenda
  theme_clean() +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 22, color = "black"),
        plot.background = element_rect(color = NA),
        legend.position = c(0.8,0.2),
        legend.text = element_text(size = 16, color = "black"))

# Criação do gráfico com barras de erro
ggplot(data_binned.TType, aes(x = mean_logRproj_rvir, y = mean_TType, color = logMgroup_bin)) +
  geom_line(aes(group = logMgroup_bin), size = 0.8) +  # Linhas conectando os pontos
  geom_point(shape = 15, size = 3) +  # Pontos quadrados
  geom_errorbar(aes(ymin = mean_TType - error_TType, ymax = mean_TType + error_TType),
                width = 0.05, size = 0.5) +  # Barras de erro
  scale_color_manual(values = c("red", "orange", "blue"),  
                     labels = c("12.3 < log Mh < 13.5",
                                "13.5 < log Mh < 14.0",
                                "14.0 < log Mh < 15.5")) +
  scale_y_continuous(breaks = pretty_breaks(n = 6)) + 
  scale_x_continuous(breaks = pretty_breaks(n = 6)) + 
  labs(x = expression(log(R[proj] / r[vir])),
       y = expression("Mean TType"),
       color = NULL) +  # Remove título da legenda
  theme_clean() +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 22, color = "black"),
        plot.background = element_rect(color = NA),
        legend.position = c(0.8,0.2),
        legend.text = element_text(size = 16, color = "black"))

# Criação do gráfico com barras de erro
ggplot(data_binned.velDisp_e, aes(x = mean_logRproj_rvir, y = mean_velDisp_e, color = logMgroup_bin)) +
  geom_line(aes(group = logMgroup_bin), size = 0.8) +  # Linhas conectando os pontos
  geom_point(shape = 15, size = 3) +  # Pontos quadrados
  geom_errorbar(aes(ymin = mean_velDisp_e - error_velDisp_e, ymax = mean_velDisp_e + error_velDisp_e),
                width = 0.05, size = 0.5) +  # Barras de erro
  scale_color_manual(values = c("red", "orange", "blue"),  
                     labels = c("12.3 < log Mh < 13.5",
                                "13.5 < log Mh < 14.0",
                                "14.0 < log Mh < 15.5")) +
  scale_y_continuous(breaks = pretty_breaks(n = 6)) + 
  scale_x_continuous(breaks = pretty_breaks(n = 6)) + 
  labs(x = expression(log(R[proj] / r[vir])),
       y = expression("Mean velDisp_e"),
       color = NULL) +  # Remove título da legenda
  theme_clean() +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 22, color = "black"),
        plot.background = element_rect(color = NA),
        legend.position = c(0.8,0.8),
        legend.text = element_text(size = 16, color = "black"))

# Modelos
fit_fSFG <- glm(SF ~ logvelDisp_e + logRproj_rvir + logMgroup, family = binomial(link = "logit"), data = data.s_cut)
fit_fLTG <- glm(LT ~ logvelDisp_e + logRproj_rvir + logMgroup, family = binomial(link = "logit"), data = data.s_cut)

summary(fit_fSFG)
summary(fit_fLTG)
