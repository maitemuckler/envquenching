setwd("~/Work/Research/")

library(combinat)
library(dplyr)
library(data.table)

input_file <- "~/Work/Research/Astronomy/Data/InputModel/inputdata_zmax0.1_Rlim2.5_Ma12.3_flag_good==1_MANGLE_logMstar_min10.5.csv"

df <- fread(input_file)
df <- subset(df, df$type == "Satellite")
#df <- subset(df, df$AGN == "AGN")
df <- df %>%
  rename(SF = SF_GSWLC)

variaveis <- c("SF",
               #"AGN",
               "logvelDisp_e",
               "logMstar",
               "B_T_r",
               "TType",
               "conc",
               "logMgroup",
               "logRproj_rvir",
               "logdens_n_Rproj0.5",
               "logdens_n_Rproj1",
               "logdens_proj_Neq1",
               "logdens_proj_Neq3",
               "logdens_proj_Neq5")

df <- df %>%
  select(all_of(variaveis))

df$SF          <- ifelse(df$SF == "Star-forming", 1, 0)
df$SF          <- as.factor(df$SF)

df$AGN         <- ifelse(df$AGN == "AGN", 1, 0)
df$AGN         <- as.factor(df$AGN)

sapply(df, class)

# Removendo SF do vetor "variáveis"
variaveis <- variaveis[-which(variaveis == "SF")]

# Listas para armazenar os valores de AIC e taxas de acerto
valores_aic  <- list()
taxas_acerto <- list()

# Gerar todas as combinações possíveis de variáveis
for (i in 1:length(variaveis)) {
  
  print(paste0(i,"/",length(variaveis)))
  
  combinacoes <- combn(variaveis, i, simplify = FALSE)
  
  # Iterar sobre cada combinação e ajustar o modelo logístico
  for (comb in combinacoes) {
    formula <- as.formula(paste("SF ~", paste(comb, collapse = " + ")))
    modelo <- glm(formula, data = df, family = binomial) # Ajustar o modelo logístico usando glm
    
    # Criar um identificador de modelo baseado nas variáveis usadas
    nome_modelo <- paste(comb, collapse = " + ")
    
    # Calcular o AIC do modelo e armazenar na lista
    valores_aic[[nome_modelo]] <- AIC(modelo)
    
    # Gerar valores preditos
    prob_predita   <- predict(modelo, type = "response")
    optimal        <- optimalCutoff(df$SF, prob_predita)[1]
    classe_predita <- ifelse(prob_predita > optimal, 1, 0)
    
    # Calcular a taxa de acerto
    acertos <- sum(classe_predita == df$SF)
    taxa_acerto <- acertos / nrow(df)
    taxas_acerto[[nome_modelo]] <- taxa_acerto
  }
}

# Exibir os valores de AIC de cada modelo com identificadores
print(valores_aic)

# Exibir as taxas de acerto de cada modelo com identificadores
print(taxas_acerto)

min_aic <- min(unlist(valores_aic))  # O menor valor de AIC
min_aic_model <- names(valores_aic)[which.min(unlist(valores_aic))]  # Modelo correspondente

min_aic
min_aic_model

max_taxa_acerto       <- max(unlist(taxas_acerto))  # A maior taxa de acerto
max_taxa_acerto_model <- names(taxas_acerto)[which.max(unlist(taxas_acerto))]  # Modelo correspondente

taxas_df <- stack(taxas_acerto)
taxas_df <- taxas_df %>%
  rename(taxas_acerto = values)

aic_df <- stack(valores_aic)
aic_df <- aic_df %>%
  rename(aic = values)

modelos <- merge(taxas_df, aic_df, by = "ind")
write.csv(modelos, "bestSubset_2.csv", row.names = F)

# ----------------------------------------
final_model <- glm(SF ~ logvelDisp_e + TType + d4000_n, 
                   family = binomial(link = "logit"), data = df)

summary(final_model)

prob_predita <- predict(final_model, type = "response")
classe_predita <- ifelse(prob_predita > 0.5, 1, 0)

# Calcular a taxa de acerto
acertos <- sum(classe_predita == df$SF)
taxa_acerto <- acertos / nrow(df)
taxa_acerto
