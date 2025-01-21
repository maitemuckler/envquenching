library(dplyr)
library(data.table)

### Definir diretórios ----
wdmain       <- "~/Work/Research/Astronomy/"
wdproject    <- paste0(wdmain, "Projects/environmental-quenching/")
wdcode       <- paste0(wdproject, "Scripts/")
wddata       <- paste0(wdmain, "Data/Raw data/SDSS/")

### Ler arquivo de entrada ----
file <- "SDSS_DR18_Legacy_MGS_QSO.csv"
df   <- fread(paste0(wddata, file))

### Ler funções ----
source(paste0(wdcode, "MyFunctions/DM_calc.R"))

### Calcular módulo da distância para os redshifts ----
DM        <- DM_calc(redshifts = df$z)
df$DM_mag <- DM

# Criar magnitudes absolutas ----
df$absPetro_u <- df$magPetro_u_kcorr - df$DM_mag
df$absPetro_g <- df$magPetro_g_kcorr - df$DM_mag
df$absPetro_r <- df$magPetro_r_kcorr - df$DM_mag
df$absPetro_i <- df$magPetro_i_kcorr - df$DM_mag
df$absPetro_z <- df$magPetro_z_kcorr - df$DM_mag

df$igal <- 1:nrow(df)

colOrder <- c(colnames(df)[73], colnames(df)[1:34], colnames(df)[67:72], colnames(df)[35:66])
df       <- df %>% select(all_of(colOrder))

write.csv(df, paste0(wddata, file), row.names = F)
