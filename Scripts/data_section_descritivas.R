library(corrplot)
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)

setwd("~/Work/Research/")

input_file <- "~/Work/Research/Astronomy/Data/InputModel/inputdata_zmax0.1_Rlim2.5_Ma12.3_flag_good==1_MANGLE_logMstar_min10.5.csv"

df <- fread(input_file)

df <- df %>%
  rename(SF = SF_GSWLC)

df <- subset(df, df$type == "Satellite")

# Variáveis que quero trabalhar
variaveis <- c("SF",
               "logvelDisp_e",
               "logMstar",
               "B_T_r",
               "TType",
               "conc",
               "d4000_n",
               "logMgroup",
               "logRproj_rvir",
               "logdens_n_Rproj0.5",
               "logdens_n_Rproj1",
               "logdens_proj_Neq1",
               "logdens_proj_Neq3",
               "logdens_proj_Neq5")

var_int <- c("SF",
             "logvelDisp_e",
             "logMstar",
             "B_T_r",
             "TType",
             "conc",
             "d4000_n")

var_ext <- c("SF",
             "logMgroup",
             "logRproj_rvir",
             "logdens_n_Rproj0.5",
             "logdens_n_Rproj1",
             "logdens_proj_Neq1",
             "logdens_proj_Neq3",
             "logdens_proj_Neq5")

df <- df %>%
  select(all_of(variaveis))

df_int <- df %>%
  select(all_of(var_int))
df_int_long <- pivot_longer(df_int, cols = logvelDisp_e:d4000_n, names_to = "Variable", values_to = "Value")

df_ext <- df %>%
  select(all_of(var_ext))
df_ext_long <- pivot_longer(df_ext, cols = logMgroup:logdens_proj_Neq5, names_to = "Variable", values_to = "Value")

# Relação com a resposta
# Variáveis intrínsecas
outliers_d4000_n  <- boxplot.stats(df_int$d4000_n)$out
outliers_conc     <- boxplot.stats(df_int$conc)$out
outliers_logMstar <- boxplot.stats(df_int$logMstar)$out

df_int_long <- df_int_long %>%
  filter(!(Variable == "d4000_n" & Value %in% outliers_d4000_n)) %>%
  filter(!(Variable == "conc" & Value %in% outliers_conc)) %>%
  filter(!(Variable == "logMstar" & Value %in% outliers_logMstar))

ggplot(df_int_long, aes(x = Value, fill = SF)) + 
  geom_density(alpha = 0.3) + 
  facet_wrap(.~Variable, scales = "free")

# Variáveis extrínsecas
ggplot(df_ext_long, aes(x = Value, fill = SF)) + 
  geom_density(alpha = 0.3) + 
  facet_wrap(.~Variable, scales = "free")

# Correlações:
SF <- df_int %>%
  filter(SF == "Star-forming")
SF_corr <- SF %>% 
  select(-SF) 
SF_corr <- cor(SF_corr)
corrplot(SF_corr, type = 'lower', diag = F,
         tl.col = "black", bg = "White", 
         addCoef.col = "black")

Q <- df_int %>%
  filter(SF == "Quiescent")
Q_corr <- Q %>% 
  select(-SF) 
Q_corr <- cor(Q_corr)
corrplot(Q_corr, type = 'lower', diag = F,
         tl.col = "black", bg = "White", 
         addCoef.col = "black")

SF <- df_ext %>%
  filter(SF == "Star-forming")
SF_corr <- SF %>% 
  select(-SF) 
SF_corr <- cor(SF_corr)
corrplot(SF_corr, type = 'lower', diag = F,
         tl.col = "black", bg = "White", 
         addCoef.col = "black")

Q <- df_ext %>%
  filter(SF == "Quiescent")
Q_corr <- Q %>% 
  select(-SF) 
Q_corr <- cor(Q_corr)
corrplot(Q_corr, type = 'lower', diag = F,
         tl.col = "black", bg = "White", 
         addCoef.col = "black")
