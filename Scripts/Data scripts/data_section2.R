## Lendo pacotes ----
library(data.table)
library(dplyr)

## Diretórios ----
wdcode <- "Scripts/"
wddata <- "/home/muckler/Work/Research/Astronomy/Data/"
wdfigs <- "Figures/"

## Definir qual tabela assignment ----
zmax    <- 0.1
Rlim    <- 2.5
Ma      <- 12.3

## Definindo input e output files ----
input_clean_file  <- "clean_letter1_sample.csv"
input_assing_file <- paste0("assignment2groups_zmax",zmax,"_Rlim",Rlim,"_Ma",Ma,".csv")
output_file <- paste0("assign_zmax", zmax, "_Ma", Ma, "_", input_clean_file)

## Lendo os dados ----
df_clean  <- fread(paste0(wddata, input_clean_file)) # 255.727 (328020 - 255727 = 72293 excluídas)
df_assign <- fread(paste0(wddata, "Assignment2groups/", input_assing_file))

# df_assign para zmax = 0.03 e Ma = 12.3: 26.937
# df_assign para zmax = 0.10 e Ma = 12.3: 

## Unir tabela clean e assignment ----
df <- merge(df_assign, df_clean, by = 'igal') 

# df merged para zmax = 0.03 e Ma = 12.3: 14.718
# df merged para zmax = 0.10 e Ma = 12.3: 

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
                            "groupID")
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

## Salvar tabela ----
write.csv(df, paste0(wddata, "Assignment2groups/", output_file), row.names = F)

