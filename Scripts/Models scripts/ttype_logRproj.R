setwd("~/Work/Research/Astronomy/Projects/environmental-quenching/Scripts/Models scripts/")

library(data.table)
library(binom)
library(ggplot2)
library(InformationValue)
library(caret)
library(dplyr)

source("~/Work/Research/Astronomy/Projects/environmental-quenching/Scripts/Themes/my_theme.R")

zmax       <- 0.03
input_data <- paste0("inputdata_zmax",zmax,"_Rlim2.5_Ma12.3_flag_good==1_MANGLE_logMstar_min10.5.csv")

df    <- fread(paste0("~/Work/Research/Astronomy/Data/EnvQuenching/inputModel/GSWLC/", input_data))
df    <- df[-which(df$logvelDisp_e < log10(50)),]
df$SF <- ifelse(df$SF_GSWLC == "Star-forming", 1, 0)
df$SF <- as.factor(df$SF)
df    <- subset(df, df$type == "Satellite")

sigma_breaks <- c(min(df$velDisp_e), 130, 150, max(df$velDisp_e)) 

# Dividir os dados em bins e calcular a média em cada bin
data_binned <- df %>%
  mutate(TType_bin = cut(TType, breaks = 50),
         logRproj_rvir_bin = cut(logRproj_rvir, breaks = 50),
         sigma_bin = cut(velDisp_e, breaks = sigma_breaks, include.lowest = TRUE)) %>%
  group_by(TType_bin, logRproj_rvir_bin, sigma_bin, SF) %>%
  summarise(median_TType = median(TType),
            median_logRproj_rvir = median(logRproj_rvir)) %>%
  ungroup()

# Plotar o gráfico
ggplot(data_binned, aes(x = median_logRproj_rvir, y = median_TType, color = factor(SF))) +
  geom_smooth(aes(group = SF)) +
  labs(x = "logRproj", y = "t-type", color = "Class") +
  scale_color_manual(values = c("red", "blue")) 

# Plotar o gráfico
ggplot(data_binned, aes(x = median_logRproj_rvir, y = median_TType, color = factor(SF))) +
  geom_smooth(aes(group = SF)) +
  facet_grid(. ~ sigma_bin) + 
  labs(x = "logRproj", y = "t-type", color = "Class") +
  scale_color_manual(values = c("red", "blue")) 
