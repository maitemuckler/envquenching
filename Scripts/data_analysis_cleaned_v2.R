
# =============================
# CONFIGURAÇÃO INICIAL
# =============================

## Pacotes necessários ----
library(data.table)
library(ggplot2)
library(scales)
library(patchwork)
library(dplyr)
library(extrafont)
library(binom)
library(InformationValue)
library(caret)
library(viridis)
library(ggthemes)

# Carregar temas e funções personalizadas
wdcode <- "Scripts/"
source(paste0(wdcode, "Themes/my_theme.R"))
source(paste0(wdcode, "MyFunctions/sigma_gal_correction.R"))
source(paste0(wdcode, "MyFunctions/SF_Q_class.R"))
source(paste0(wdcode, "MyFunctions/distance_to_line_calc.R"))

# Configurações gerais
set.seed(123)
figs_width  <- 14
figs_height <- 8

# Diretórios principais
wddata <- "/home/muckler/Work/Research/Astronomy/Data/"
wdfigs <- "Figures/"

# =============================
# PARÂMETROS DO ESTUDO
# =============================

zmax <- 0.1
Rlim <- 2.5
Ma   <- 12.3
catalog <- "GSWLC"

SF_name <- ifelse(catalog == "GSWLC", "SF_GSWLC", "SF_MPAJHU")
mass_complete <- "y"

input_file_clean   <- "letter1_sample.csv"
input_file_assign  <- paste0("assignment2groups_zmax", zmax, "_Rlim", Rlim, "_Ma", Ma, ".csv")
output_file_clean  <- paste0("clean_", input_file_clean)
output_file_assign <- paste0("assign_zmax", zmax, "_Ma", Ma, "_", input_file_clean)

# =============================
# LEITURA E PRÉ-PROCESSAMENTO DOS DADOS
# =============================

## Leitura dos dados ----
df_clean  <- fread(paste0(wddata, input_file_clean))
df_assign <- fread(paste0(wddata, "Assignment2groups/", input_file_assign))

## Mesclagem das tabelas ----
df <- merge(df_assign, df_clean, by = 'igal')

# Removendo colunas duplicadas do merge
colunas_para_remover <- grep("\.y$", names(df), value = TRUE)
df <- df %>% select(-all_of(colunas_para_remover))
rm(colunas_para_remover)

# Renomeação de colunas
colunas_para_renomear <- grep("\.x$", names(df), value = TRUE)
colnames(df)[which(colnames(df) %in% colunas_para_renomear)] <- gsub("\.x$", "", colunas_para_renomear)
rm(colunas_para_renomear)

## Criando novas variáveis ----
df$logRproj_rvir <- log10(df$Rproj_rvir)
df$type <- factor(ifelse(df$Rproj_rvir > 0, "Satellite", "Central"))

# Convertendo variáveis para fatores
fatores <- c("Class", "bptClass", "flag_good", "flag_sed", "flag_uv", "flag_midir", 
             "type", "SF_GSWLC", "SF_MPAJHU", "TType_label", "groupID")
df <- df %>% mutate(across(all_of(fatores), factor))
rm(fatores)

## Filtragem ----
df <- df %>% filter(flag_good == 1)

## Correção do MANGLE ----
mangle_L  <- fread("~/Work/Research/Astronomy/Data/Mangle/groups_fmangle_L.csv")
groups_rm <- mangle_L$groupID_L[which(mangle_L$f_mangle_20rvir_L < 0.9)]
df <- subset(df, !(df$groupID %in% groups_rm))
rm(mangle_L)

## Completeza em massa estelar ----
zlim_logMstar <- fread("~/Work/Research/Astronomy/Data/zlim_logMstar/logM_zlim_MAITE.csv")
logMstar_min  <- zlim_logMstar$logM[which(zlim_logMstar$z_max_seq == zmax)]
df_cutmass    <- subset(df, df$lgm_tot_p50 >= logMstar_min)
rm(zlim_logMstar)

## Salvamento dos arquivos filtrados ----
write.csv(df, paste0(wddata, "inputModel/", catalog, "/input_filtered.csv"), row.names = F)
write.csv(df_cutmass, paste0(wddata, "inputModel/", catalog, "/input_filtered_mass.csv"), row.names = F)

# =============================
# MODELAGEM ESTATÍSTICA
# =============================

## Ajuste de modelos logísticos ----

if(mass_complete == "y") {
  fit_fSFG <- glm(SF ~ logvelDisp_e + logRproj_rvir, family = binomial(link = "logit"), data = df)
  fit_fLTG <- glm(LT ~ logvelDisp_e + logRproj_rvir, family = binomial(link = "logit"), data = df)
} else {
  df$Vmax[df$Vmax == 0] <- 1e-6
  df$weights <- 1 / df$Vmax
  df$weights <- df$weights / max(df$weights)
  
  fit_fSFG <- glm(SF ~ logvelDisp_e + logRproj_rvir, family = binomial(link = "logit"), data = df, weights = weights)
  fit_fLTG <- glm(LT ~ logvelDisp_e + logRproj_rvir, family = binomial(link = "logit"), data = df, weights = weights)
}

summary(fit_fSFG)
summary(fit_fLTG)

# =============================
# VISUALIZAÇÃO
# =============================

ggplot(df, aes(x = logvelDisp_e, y = logRproj_rvir, color = factor(SF))) +
  geom_point(alpha = 0.5) +
  labs(x = "Log VelDisp", y = "Log Rproj", color = "Star-Forming") +
  theme_minimal()

ggsave(paste0(wdfigs, "model_visualization.png"), width = figs_width, height = figs_height, dpi = 300)
