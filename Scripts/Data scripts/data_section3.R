## Lendo pacotes ----
library(data.table)

## Diretórios ----
wdcode <- "Scripts/"
wddata <- "/home/muckler/Work/Research/Astronomy/Data/"
wdfigs <- "Figures/"

## Definir qual tabela assignment ----
zmax    <- 0.03
Rlim    <- 2.5
Ma      <- 12.3
catalog <- "GSWLC"

## Definindo input e output files ----
input_file          <- paste0("assign_zmax", zmax, "_Ma", Ma, "_clean_letter1_sample.csv")
output_file         <- paste0("inputdata_zmax", zmax, "_Rlim", Rlim, "_Ma", Ma, "_flag_good==1_MANGLE.csv")
output_file_cutMass <- paste0("inputdata_zmax", zmax, "_Rlim", Rlim, "_Ma", Ma, "_flag_good==1_MANGLE_logMstar_min10.5.csv")

## Lendo os dados ----
df  <- fread(paste0(wddata, "Assignment2groups/", input_file)) 

## Corrigindo MANGLE ----
mangle_L  <- fread("~/Work/Research/Astronomy/Data/Mangle/groups_fmangle_L.csv")
groups_rm <- mangle_L$groupID_L[which(mangle_L$f_mangle_20rvir_L < 0.9)]
'%ni%'    <- Negate("%in%")
df <- subset(df, df$groupID %ni% groups_rm) # 181.249
rm(mangle_L)

## Filtrando para flag_good == 1 ----
df <- df %>%
  filter(flag_good == 1) # 170.926

## Completeza em massa estelar ----
zlim_logMstar <- fread("~/Work/Research/Astronomy/Data/zlim_logMstar/logM_zlim_MAITE.csv")
logMstar_min  <- zlim_logMstar$logM[which(zlim_logMstar$z_max_seq == zmax)]
ngals_bin     <- zlim_logMstar$Ngals_bin[which(zlim_logMstar$z_max_seq == zmax)]

#df_cutmass    <- subset(df, df$logMstar >= logMstar_min) # 75.140
df_cutmass    <- subset(df, df$lgm_tot_p50 >= logMstar_min) 
rm(zlim_logMstar)

## Verificando valores infinitos e NA ----
names(df_cutmass)[sapply(df_cutmass, function(x) any(is.na(x)))]
names(df_cutmass)[sapply(df_cutmass, function(x) any(is.infinite(x)))]

## Salvando tabelas ----
write.csv(df, paste0(wddata, "inputModel/", catalog, "/", output_file), row.names = F)
write.csv(df_cutmass, paste0(wddata, "inputModel/", catalog, "/", output_file_cutMass) , row.names = F)