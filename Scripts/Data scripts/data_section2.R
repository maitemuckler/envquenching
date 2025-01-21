## Lendo pacotes ----
library(data.table)
library(dplyr)

## Diretórios ----
wdcode <- "~/Work/Research/Astronomy/Projects/environmental-quenching/Scripts/"
wddata <- "~/Work/Research/Astronomy/Data/environmental-quenching-data/"
wdfigs <- "~/Work/Research/Astronomy/Projects/environmental-quenching/Figures/"

## Definir qual tabela assignment ----
zmax    <- 0.07
Rlim    <- 2.5
Ma      <- 14

## Definindo input e output files ----
input_clean_file  <- "clean_SDSS_DR18_Legacy_MGS_QSO_localDensity+GSWLC+Simard11+DS18.csv" 
input_assing_file <- paste0("assignment2groups_zmax",zmax,"_Rlim",Rlim,"_Ma",Ma,".csv")

output_file <- paste0("assign_zmax", zmax, "_Ma", Ma, "_", input_clean_file)

## Lendo os dados ----
df_clean  <- fread(paste0(wddata, input_clean_file)) # 255.688 (328.020 - 255.688 = 72.332 excluídas)
df_assign <- fread(paste0(wddata, "Assignment2groups/", input_assing_file))

# Para zmax = 0.1 e Ma = 12.3:   297.555 (320.020 - 297.55 = 30.465 galáxias que não estão em grupos)

## Unir tabela clean e assignment ----
df <- merge(df_assign, df_clean, by = 'igal') 

# Para zmax = 0.1 e Ma = 12.3: 232.268

#  A alegação é:
# "Das 255.688 galáxias da nossa amostra, 23.420 não foram atribuidas a nenhum grupo."
#  Questão: São galáxias isoladas ou apenas galáxias que não estão em nenhum grupo do Lim?

## Remover essas colunas repetidas no merge ----
colunas_para_remover <- grep("\\.y$", names(df), value = TRUE)

df <- df %>% 
  select(-all_of(colunas_para_remover))

rm(colunas_para_remover)

## Renomear colunas repetidas no merge ----
colunas_para_renomear <- grep("\\.x$", names(df), value = TRUE)

colnames(df)[which(colnames(df) %in% colunas_para_renomear)] <- gsub("\\.x$", "", colunas_para_renomear)

rm(colunas_para_renomear)

## Criar variáveis (pós-assignment) ----

# log10 da distância clustercentrica:
df$logRproj_rvir <- log10(df$Rproj_rvir) # log(Rproj_rvir)

# Indicadora se é satélite ou central:
df$type <- ifelse(df$Rproj_rvir > 0, "Satellite", "Central")
df$type <- as.factor(df$type)

# log10(Mstar/Mhalo): 
df$Mstar_Mhalo <- df$logMstar - df$logMgroup

## Criar colunas para variáveis a respeito das centrais ----

df$groupID <- as.factor(df$groupID)

# Colunas que quero salvar a informação da central 
colsCentral <- c("groupID", 
                 "logMstar", 
                 "logSFR_SED", 
                 "logvelDisp_e", 
                 "vlos_vvir", 
                 "conc", 
                 "d4000_n",
                 "logSigma_SFR", 
                 "distLine_GSWLC", 
                 "SF_GSWLC", 
                 "AGN",
                 "P_disk", 
                 "P_edge_on", 
                 "P_bar_GZ2", 
                 "P_bar_Nair10", 
                 "P_merg", 
                 "P_bulge", 
                 "B_T_r", 
                 "e", 
                 "P_cigar", 
                 "TType", 
                 "P_S0", 
                 "TType_label")

colsCentral_names <- colsCentral
colsCentral_names[-1] <- paste0(colsCentral[-1], "_central")

aux <- df %>%
  filter(type == "Central") %>%
  select(all_of(colsCentral))

colnames(aux) <- colsCentral_names

df <- merge(df, aux, by = "groupID", all.x = T)

rowIndex <- which(colnames(df) %in% colsCentral_names)
rowIndex <- rowIndex[-1]
nrow(df[rowSums(is.na(df[,..rowIndex])) == ncol(df[,..rowIndex]), ..rowIndex]) # Sem informação de central

rm(aux, df_assign, df_clean)

## Converter variáveis para classes corretas ----

lapply(df, class) # Verificar o tipo das variáveis

# Variáveis que quero transformar para factor:
colunas_para_converter <- c("Class",
                            "bptClass",
                            "flag_good",
                            "flag_sed",
                            "flag_uv",
                            "flag_midir",
                            "type",
                            "SF_GSWLC",
                            "SF_MPAJHU",
                            "TType_label",
                            "AGN",
                            "groupID", 
                            "SF_GSWLC_central", 
                            "AGN_central",
                            "TType_label_central")
df <- df %>%
  mutate(across(all_of(colunas_para_converter), factor))

rm(colunas_para_converter)

## Adicionar label em algumas variáveis do tipo factor ---

# ---- bptClass:

# Emission line classification based on the BPT diagram using the methodology described in Brinchmann et al (2004). 
# -1 means unclassifiable, 
# 1 is star-forming, 
# 2 means low S/N star-forming, 
# 3 is composite, 
# 4 AGN (excluding liners) and 
# 5 is a low S/N LINER.

levels(df$bptClass) <- c("-1" = "unclassifiable", 
                         "1"  = "star-forming",
                         "2"  = "low S/N star-forming",
                         "3"  = "composite",
                         "4"  = "AGN (excluding liners)",
                         "5"  = "low S/N LINER")

# ---- flag_sed:

# SED fitting flag 
# 0 = OK, 
# 1 = broad-line spectrum, 
# 2 = X²_r > 30, (reduced chi-square) 
# 5 = missing SDSS photometry

levels(df$flag_sed) <- c("0" = "OK",
                         "1" = "broad-line spectrum",
                         "2" = "X²_r > 30",
                         "3" = "Não encontrei label",
                         "5" = "missing SDSS photometry")

# ---- flag_uv:

# UV (GALEX) flag 
# 0 = no UV; 
# 1 = FUV only; 
# 2 = NUV only; 
# 3 = both

levels(df$flag_uv) <- c("0" = "no UV",
                        "1" = "FUV only",
                        "2" = "NUV only",
                        "3" = "both")

# ---- flag_midir:

# Mid-IR (unWISE) flag 
# 0 = no mid-IR, 
# 1 = LIR based on 12 μm, 
# 2 = LIR based on 22 μm;
# 5 = LIR corrected for mid-IR AGN emission

levels(df$flag_midir) <- c("0" = "no mid-IR",
                           "1" = "LIR based on 12 μm",
                           "2" = "LIR based on 22 μm",
                           "5" = "LIR corrected for mid-IR AGN emission",
                           "6" = "Não encontrei label 6",
                           "7" = "Não encontrei label 7")

## Verificando valores NA e inf ----

names(df)[sapply(df, function(x) any(is.na(x)))]
names(df)[sapply(df, function(x) any(is.infinite(x)))]

# logRproj_rvir        <- infinito porque era central, então Rproj = 0 -> logRproj = Inf.
# logO3Hb e logN2Ha    <- alguma divisão por zero.
# logvelDisp_e_central <- a central deve ter velDisp_e = 0 -> logvelDisp_e = Inf.
# Para logvelDisp_e_central vou susbtituir o Inf por NA.

# Central
df$logvelDisp_e_central[which(is.infinite(df$logvelDisp_e_central))] <- NA

## Salvar tabela ----
write.csv(df, paste0(wddata, "assign/", output_file), row.names = F)
