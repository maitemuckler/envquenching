## Pacotes ----
library(data.table)
library(dplyr)

## Diretórios ----
wdcode <- "Scripts/"
wddata <- "/home/muckler/Work/Research/Astronomy/Data/"
wdfigs <- "Figures/"

## Definição de parâmetros ----
zmax <- 0.03
Rlim <- 2.5
Ma   <- 12.3

## Definição de arquivos ----
input_clean_file  <- "clean_letter1_sample.csv"
input_assing_file <- paste0("assignment2groups_zmax",zmax,"_Rlim",Rlim,"_Ma",Ma,".csv")
output_file <- paste0("assign_zmax", zmax, "_Ma", Ma, "_", input_clean_file)

## Leitura de dados ----
df_clean  <- fread(paste0(wddata, input_clean_file))
df_assign <- fread(paste0(wddata, "Assignment2groups/", input_assing_file))

## Unindo tabelas ----
df <- merge(df_assign, df_clean, by = "igal")

## Removendo colunas repetidas ----
colunas_para_remover <- grep("\\.y$", names(df), value = TRUE)
df <- df %>% select(-all_of(colunas_para_remover))

## Renomeando colunas repetidas ----
colunas_para_renomear <- grep("\\.x$", names(df), value = TRUE)
colnames(df)[which(colnames(df) %in% colunas_para_renomear)] <- gsub("\\.x$", "", colunas_para_renomear)

## Criando variáveis pós-assignment ----
df$logRproj_rvir <- log10(df$Rproj_rvir)
df$type <- factor(ifelse(df$Rproj_rvir > 0, "Satellite", "Central"))

## Conversão de variáveis ----
lapply(df, class)
colunas_para_converter <- c("Class", "bptClass", "flag_good", "flag_sed", "flag_uv", "flag_midir", "type", "SF_GSWLC", "SF_MPAJHU", "TType_label", "morph_char", "groupID")
df <- df %>% mutate(across(all_of(colunas_para_converter), factor))

## Adicionando labels para variáveis factor ----
df$bptClass <- factor(df$bptClass, 
                      labels = c("unclassifiable", "star-forming", "low S/N star-forming", "composite", "AGN (excluding liners)", "low S/N LINER"))
df$flag_sed <- factor(df$flag_sed, 
                      levels = c(0, 1, 2, 3, 5), 
                      labels = c("OK", "broad-line spectrum", "X²_r > 30", "Não encontrei label", "missing SDSS photometry"))
df$flag_uv  <- factor(df$flag_uv, 
                      labels = c("no UV", "FUV only", "NUV only", "both"))
df$flag_midir <- factor(df$flag_midir, 
                        labels = c("no mid-IR", "LIR based on 12 μm", "LIR based on 22 μm", "LIR corrected for mid-IR AGN emission", "Não encontrei label 6", "Não encontrei label 7"))

## Verificação de NAs e valores infinitos ----
names(df)[sapply(df, function(x) any(is.na(x)))]
names(df)[sapply(df, function(x) any(is.infinite(x)))]

## Salvando tabela ----
write.csv(df, paste0(wddata, "Assignment2groups/", output_file), row.names = F)
