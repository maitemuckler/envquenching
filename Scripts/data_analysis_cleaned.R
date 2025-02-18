
## Pacotes necessários ----
library(data.table)
library(ggplot2)
library(scales)
library(patchwork)
library(dplyr)
library(extrafont)

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

## Parâmetros ----
zmax    <- 0.1
Rlim    <- 2.5
Ma      <- 12.3
catalog <- "GSWLC"

## Definição de arquivos ----
input_file_clean   <- "letter1_sample.csv"
input_file_assign  <- paste0("assignment2groups_zmax", zmax, "_Rlim", Rlim, "_Ma", Ma, ".csv")
output_file_clean  <- paste0("clean_", input_file_clean)
output_file_assign <- paste0("assign_zmax", zmax, "_Ma", Ma, "_", input_file_clean)
output_file_final  <- paste0("inputdata_zmax", zmax, "_Rlim", Rlim, "_Ma", Ma, "_flag_good==1_MANGLE.csv")
output_file_mass   <- paste0("inputdata_zmax", zmax, "_Rlim", Rlim, "_Ma", Ma, "_flag_good==1_MANGLE_logMstar_min10.5.csv")

## Leitura dos dados ----
df_clean  <- fread(paste0(wddata, input_file_clean))
df_assign <- fread(paste0(wddata, "Assignment2groups/", input_file_assign))

## Unificação das tabelas ----
df <- merge(df_assign, df_clean, by = 'igal')

## Remoção de colunas duplicadas do merge ----
colunas_para_remover <- grep("\.y$", names(df), value = TRUE)
df <- df %>% select(-all_of(colunas_para_remover))
rm(colunas_para_remover)

## Renomeação de colunas duplicadas ----
colunas_para_renomear <- grep("\.x$", names(df), value = TRUE)
colnames(df)[which(colnames(df) %in% colunas_para_renomear)] <- gsub("\.x$", "", colunas_para_renomear)
rm(colunas_para_renomear)

## Criação de variáveis adicionais ----
df$logRproj_rvir <- log10(df$Rproj_rvir)
df$type <- factor(ifelse(df$Rproj_rvir > 0, "Satellite", "Central"))

## Conversão de variáveis para fatores ----
fatores <- c("Class", "bptClass", "flag_good", "flag_sed", "flag_uv", "flag_midir", "type", 
             "SF_GSWLC", "SF_MPAJHU", "TType_label", "groupID")
df <- df %>% mutate(across(all_of(fatores), factor))
rm(fatores)

## Filtragem para flag_good == 1 ----
df <- df %>% filter(flag_good == 1)

## Remoção de grupos problemáticos do MANGLE ----
mangle_L  <- fread("~/Work/Research/Astronomy/Data/Mangle/groups_fmangle_L.csv")
groups_rm <- mangle_L$groupID_L[which(mangle_L$f_mangle_20rvir_L < 0.9)]
df <- subset(df, !(df$groupID %in% groups_rm))
rm(mangle_L)

## Completeza em massa estelar ----
zlim_logMstar <- fread("~/Work/Research/Astronomy/Data/zlim_logMstar/logM_zlim_MAITE.csv")
logMstar_min  <- zlim_logMstar$logM[which(zlim_logMstar$z_max_seq == zmax)]
df_cutmass    <- subset(df, df$lgm_tot_p50 >= logMstar_min)
rm(zlim_logMstar)

## Salvamento dos arquivos finais ----
write.csv(df, paste0(wddata, "inputModel/", catalog, "/", output_file_final), row.names = F)
write.csv(df_cutmass, paste0(wddata, "inputModel/", catalog, "/", output_file_mass), row.names = F)
