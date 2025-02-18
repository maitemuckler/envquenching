setwd("~/Work/Research/Astronomy/Projects/envquenching/")

library(data.table)
library(binom)
library(ggplot2)
library(InformationValue)
library(caret)
library(dplyr)
library(scales)
library(viridis)
library(ggthemes) # theme_clean

source("Scripts/Themes/my_theme.R")

wddata <- "~/Work/Research/Astronomy/Data/"

zmax          <- 0.03
catalog       <- "GSWLC"
TType_lim     <- 0
critval       <- 1.96 
Ma            <- 12.3

SF_name    <- ifelse(catalog == "GSWLC", "SF_GSWLC", "SF_MPAJHU")

input_data <- paste0("/inputdata_zmax", zmax, "_Rlim2.5_Ma", Ma, "_flag_good==1_MANGLE_logMstar_min10.5.csv")

df    <- fread(paste0(wddata, "inputModel/", catalog, input_data))

colnames(df)[which(colnames(df) == SF_name)] <- "SF_char"

df$SF_char  <- as.factor(df$SF_char)

df$SF          <- ifelse(df$SF_char == "Star-forming", 1, 0)

df$LT <- ifelse(df$TType >= TType_lim, 1, 0)
df$LT <- as.factor(df$LT)

data.s <- subset(df, df$type == "Satellite")
data.c <- subset(df, df$type == "Central")

#data.s <- subset(data.s, data.s$velDisp_e <= 200)

# Descritivas ------------------------------------------------------------------------------------------------

# Defina aqui os limites dos bins de sigma
sigma_breaks <- c(min(data.s$velDisp_e), 130, 150, max(data.s$velDisp_e))  
sigma_breaks <- c(min(data.s$velDisp_e), 60, 100, max(data.s$velDisp_e))  

# Dividir os dados em bins e calcular a média em cada bin
data_binned.s <- data.s %>%
  mutate(TType_bin = cut(TType, breaks = 100),
         logRproj_rvir_bin = cut(logRproj_rvir, breaks = 100),
         sigma_bin = cut(velDisp_e, breaks = sigma_breaks, include.lowest = TRUE)) %>%
  group_by(TType_bin, logRproj_rvir_bin, sigma_bin, SF) %>%
  summarise(median_TType = median(TType),
            median_logRproj_rvir = median(logRproj_rvir)) %>%
  ungroup()

# Plotar o gráfico

# Valores medianos de logRproj e T-Type por bin de sigma
ggplot(data_binned.s, aes(x = median_logRproj_rvir, y = median_TType, color = factor(SF))) +
  geom_smooth(aes(group = SF))          +
  facet_grid(. ~ sigma_bin) + 
  labs(x = "Median logRproj", y = "Median T-Type", color = "Class") +
  scale_color_manual(values = c("red", "blue")) 

# Valores medianos de logRproj e T-Type geral
ggplot(data_binned.s, aes(x = median_logRproj_rvir, y = median_TType, color = factor(SF))) +
  geom_smooth(aes(group = SF)) +
  labs(x = "Median logRproj", y = "Median T-Type", color = "Class") +
  scale_color_manual(values = c("red", "blue")) 

# Modelos ------------------------------------------------------------------------------------------------

fit_fSFG <- glm(SF ~ logvelDisp_e + logRproj_rvir, family = binomial(link = "logit"), data = data.s)
fit_fLTG <- glm(LT ~ logvelDisp_e + logRproj_rvir, family = binomial(link = "logit"), data = data.s)

summary(fit_fSFG)
summary(fit_fLTG)

importance_fSFG <- as.data.frame(varImp(fit_fSFG))
importance_fSFG <- data.frame(overall = importance_fSFG$Overall, names = rownames(importance_fSFG))
importance_fSFG[order(importance_fSFG$overall, decreasing = T),]

importance_fLTG <- as.data.frame(varImp(fit_fLTG))
importance_fLTG <- data.frame(overall = importance_fLTG$Overall, names = rownames(importance_fLTG))
importance_fLTG[order(importance_fLTG$overall, decreasing = T),]

data.s$pred_fSFG <- predict(fit_fSFG, data.s, type = "response")
data.s$pred_fLTG <- predict(fit_fLTG, data.s, type = "response")

optimal_fSFG <- optimalCutoff(data.s$SF, data.s$pred_fSFG)[1]
optimal_fLTG <- optimalCutoff(data.s$LT, data.s$pred_fLTG)[1]

misClassError(data.s$SF, data.s$pred_fSFG, threshold = optimal_fSFG)
misClassError(data.s$LT, data.s$pred_fLTG, threshold = optimal_fLTG)

# Matriz de confusão

# Converter probabilidades em classes (threshold padrão: 0.5)
data.s$pred_class_fSFG <- ifelse(data.s$pred_fSFG >= optimal_fSFG, "Star-forming", "Quiescent")  # SF = Star-Forming, Q = Quiescent
data.s$pred_class_fLTG <- ifelse(data.s$pred_fLTG >= optimal_fLTG, "Late-type", "Early-type")  # LTG = Late-Type, NT = Not Late-Type

# Gerar as matrizes de confusão para cada modelo
conf_matrix_fSFG <- confusionMatrix(factor(data.s$pred_class_fSFG), factor(data.s$SF_char))
conf_matrix_fLTG <- confusionMatrix(factor(data.s$pred_class_fLTG), factor(data.s$morph_char))

print(conf_matrix_fSFG)
print(conf_matrix_fLTG)

# Extraindo métricas de desempenho para fSFG
acc_fSFG <- conf_matrix_fSFG$overall["Accuracy"]
sensitivity_fSFG <- conf_matrix_fSFG$byClass["Sensitivity"]
specificity_fSFG <- conf_matrix_fSFG$byClass["Specificity"]

# Extraindo métricas de desempenho para fLTG
acc_fLTG <- conf_matrix_fLTG$overall["Accuracy"]
sensitivity_fLTG <- conf_matrix_fLTG$byClass["Sensitivity"]
specificity_fLTG <- conf_matrix_fLTG$byClass["Specificity"]

# Exibir resultados
cat("\nDesempenho do modelo para fSFG:\n")
cat("Acurácia:", acc_fSFG, "\nSensibilidade:", sensitivity_fSFG, "\nEspecificidade:", specificity_fSFG, "\n")

cat("\nDesempenho do modelo para fLTG:\n")
cat("Acurácia:", acc_fLTG, "\nSensibilidade:", sensitivity_fLTG, "\nEspecificidade:", specificity_fLTG, "\n")

library(pROC)

# Curva ROC para fSFG
roc_fSFG <- roc(data.s$SF, data.s$pred_fSFG)
plot(roc_fSFG, main="ROC Curve - fSFG", col="blue", lwd=2)
cat("AUC para fSFG:", auc(roc_fSFG), "\n")

# Curva ROC para fLTG
roc_fLTG <- roc(data.s$LT, data.s$pred_fLTG)
plot(roc_fLTG, main="ROC Curve - fLTG", col="red", lwd=2)
cat("AUC para fLTG:", auc(roc_fLTG), "\n")

# ODDS RATIO

# 1️⃣ Extrair coeficientes e calcular Odds Ratio para fSFG
odds_fSFG <- exp(coef(fit_fSFG))
conf_int_fSFG <- exp(confint(fit_fSFG))  # Intervalo de confiança 95%

# Criar um dataframe organizado
odds_table_fSFG <- data.frame(
  Variable = names(odds_fSFG),
  OR = odds_fSFG,
  Lower95 = conf_int_fSFG[,1],
  Upper95 = conf_int_fSFG[,2]
)

# Ordenar pelo maior efeito
odds_table_fSFG <- odds_table_fSFG[order(-odds_table_fSFG$OR),]

cat("\nOdds Ratios para fSFG:\n")
print(odds_table_fSFG)

# 2️⃣ Extrair coeficientes e calcular Odds Ratio para fLTG
odds_fLTG <- exp(coef(fit_fLTG))
conf_int_fLTG <- exp(confint(fit_fLTG))

# Criar um dataframe organizado
odds_table_fLTG <- data.frame(
  Variable = names(odds_fLTG),
  OR = odds_fLTG,
  Lower95 = conf_int_fLTG[,1],
  Upper95 = conf_int_fLTG[,2]
)

# Ordenar pelo maior efeito
odds_table_fLTG <- odds_table_fLTG[order(-odds_table_fLTG$OR),]

cat("\nOdds Ratios para fLTG:\n")
print(odds_table_fLTG)

# Teste de Multicolinearidade (VIF)
# Instalar pacote se necessário
library(car)

# Testar multicolinearidade (VIF) para cada modelo
vif_fSFG <- vif(fit_fSFG)
vif_fLTG <- vif(fit_fLTG)

cat("\nVIF para fSFG:\n")
print(vif_fSFG)

cat("\nVIF para fLTG:\n")
print(vif_fLTG)

# Ótimo! Seus valores de VIF estão todos próximos de 1, o que significa que não há multicolinearidade significativa entre as variáveis.
# Isso é um bom sinal, pois indica que os coeficientes do seu modelo não estão sendo distorcidos por correlações excessivas entre preditores.

# Verificação dos Resíduos

# Gráficos de resíduos para verificar padrões de não linearidade
par(mfrow = c(1,2))  # Criar dois gráficos lado a lado

# Resíduos vs valores ajustados para fSFG
plot(fitted(fit_fSFG), residuals(fit_fSFG), 
     main="Resíduos vs Ajustados - fSFG",
     xlab="Valores Ajustados", ylab="Resíduos", pch=16, col="blue")
abline(h=0, col="red", lty=2)

# Resíduos vs valores ajustados para fLTG
plot(fitted(fit_fLTG), residuals(fit_fLTG), 
     main="Resíduos vs Ajustados - fLTG",
     xlab="Valores Ajustados", ylab="Resíduos", pch=16, col="blue")
abline(h=0, col="red", lty=2)

#  O modelo parece bem especificado no sentido de que não há padrões aleatórios indicativos de erro grosseiro.
# No entanto, os resíduos seguem uma curva, o que pode indicar que um termo não linear ou uma transformação das variáveis 
# poderia melhorar o ajuste.

# Dado esse comportamento, aqui estão algumas análises que podemos fazer para melhorar a modelagem:
# 1️⃣ Analisar resíduos padronizados (Pearson ou Deviance)
# O que você gerou são resíduos brutos. Podemos verificar resíduos de Pearson ou Deviance residuals, que são mais apropriados 
# para modelos logísticos.
# 2️⃣ Testar a inclusão de um termo quadrático ou interação
# Se houver indícios de que as variáveis não influenciam de forma linear, podemos adicionar um termo quadrático.
# 3️⃣ Analisar a influência de pontos extremos
# Algumas observações podem estar influenciando fortemente o modelo.

# Código para Resíduos de Pearson e Deviance

par(mfrow = c(1,2))  # Dois gráficos lado a lado

# Resíduos de Pearson
plot(fitted(fit_fSFG), residuals(fit_fSFG, type = "pearson"), 
     main="Resíduos de Pearson - fSFG",
     xlab="Valores Ajustados", ylab="Resíduos", pch=16, col="blue")
abline(h=0, col="red", lty=2)

plot(fitted(fit_fLTG), residuals(fit_fLTG, type = "pearson"), 
     main="Resíduos de Pearson - fLTG",
     xlab="Valores Ajustados", ylab="Resíduos", pch=16, col="blue")
abline(h=0, col="red", lty=2)

# Resíduos Deviance
par(mfrow = c(1,2))  # Resetando o layout dos gráficos

plot(fitted(fit_fSFG), residuals(fit_fSFG, type = "deviance"), 
     main="Resíduos de Deviance - fSFG",
     xlab="Valores Ajustados", ylab="Resíduos", pch=16, col="blue")
abline(h=0, col="red", lty=2)

plot(fitted(fit_fLTG), residuals(fit_fLTG, type = "deviance"), 
     main="Resíduos de Deviance - fLTG",
     xlab="Valores Ajustados", ylab="Resíduos", pch=16, col="blue")
abline(h=0, col="red", lty=2)

# Os resíduos de Deviance ainda apresentam o mesmo padrão curvo dos resíduos brutos, 
# o que indica que o modelo logístico pode estar capturando a tendência geral, mas a relação entre as 
# variáveis preditoras e a resposta pode não ser perfeitamente linear. Esse comportamento sugere que um termo 
# não linear ou uma interação pode melhorar o ajuste do modelo.

# 1️⃣ Adicionar termo quadrático ao modelo
fit_fSFG_quad <- glm(SF ~ logvelDisp_e + I(logvelDisp_e^2) + logRproj_rvir, 
                     family = binomial(link = "logit"), data = data.s)
fit_fLTG_quad <- glm(LT ~ logvelDisp_e + I(logvelDisp_e^2) + logRproj_rvir, 
                     family = binomial(link = "logit"), data = data.s)

# 2️⃣ Adicionar interação entre variáveis
fit_fSFG_inter <- glm(SF ~ logvelDisp_e * logRproj_rvir, 
                      family = binomial(link = "logit"), data = data.s)
fit_fLTG_inter <- glm(LT ~ logvelDisp_e * logRproj_rvir, 
                      family = binomial(link = "logit"), data = data.s)

# 3️⃣ Comparar modelos usando ANOVA
anova(fit_fSFG, fit_fSFG_quad, test="Chisq")
anova(fit_fLTG, fit_fLTG_quad, test="Chisq")

anova(fit_fSFG, fit_fSFG_inter, test="Chisq")
anova(fit_fLTG, fit_fLTG_inter, test="Chisq")

# Teste de Interação Entre Variáveis (Opcional)

# Criar novos modelos com interação
fit_fSFG_inter <- glm(SF ~ logvelDisp_e * logRproj_rvir, 
                      family = binomial(link = "logit"), data = data.s)
fit_fLTG_inter <- glm(LT ~ logvelDisp_e * logRproj_rvir, 
                      family = binomial(link = "logit"), data = data.s)

# Comparar com os modelos originais usando ANOVA
anova(fit_fSFG, fit_fSFG_inter, test="Chisq")
anova(fit_fLTG, fit_fLTG_inter, test="Chisq")

if (!require(pROC)) install.packages("pROC", dependencies = TRUE)
library(pROC)

# Certificar que as colunas de resposta são numéricas (0 e 1)
data.s$SF <- as.numeric(as.character(data.s$SF))
data.s$LT <- as.numeric(as.character(data.s$LT))

# Gerar predições com o modelo logístico final
data.s$pred_fSFG_full <- predict(fit_fSFG_full, newdata = data.s, type = "response")
data.s$pred_fLTG_full <- predict(fit_fLTG_full, newdata = data.s, type = "response")

# Remover NAs, se houver
data.s <- na.omit(data.s)

# Verificar se as predições agora existem
summary(data.s$pred_fSFG_full)
summary(data.s$pred_fLTG_full)

# Criar Curva ROC para fSFG
roc_fSFG <- roc(data.s$SF, data.s$pred_fSFG_full)
plot(roc_fSFG, main="Curva ROC - fSFG", col="blue", lwd=2)
cat("AUC para fSFG:", auc(roc_fSFG), "\n")

# Criar Curva ROC para fLTG
roc_fLTG <- roc(data.s$LT, data.s$pred_fLTG_full)
plot(roc_fLTG, main="Curva ROC - fLTG", col="red", lwd=2)
cat("AUC para fLTG:", auc(roc_fLTG), "\n")

# Código para Criar o Gráfico de Odds Ratios

# Extrair coeficientes do modelo logístico e calcular Odds Ratios
odds_fSFG <- exp(coef(fit_fSFG))
conf_int_fSFG <- exp(confint(fit_fSFG))

# Criar tabela para fSFG
odds_table_fSFG <- data.frame(
  Variable = names(odds_fSFG),
  OR = odds_fSFG,
  Lower95 = conf_int_fSFG[,1],
  Upper95 = conf_int_fSFG[,2]
)

# Repetir para fLTG
odds_fLTG <- exp(coef(fit_fLTG))
conf_int_fLTG <- exp(confint(fit_fLTG))

# Criar tabela para fLTG
odds_table_fLTG <- data.frame(
  Variable = names(odds_fLTG),
  OR = odds_fLTG,
  Lower95 = conf_int_fLTG[,1],
  Upper95 = conf_int_fLTG[,2]
)

# Adicionar rótulos para identificar os modelos
odds_table_fSFG$Model <- "fSFG"
odds_table_fLTG$Model <- "fLTG"

# Unir os dois dataframes
odds_table_all <- rbind(odds_table_fSFG, odds_table_fLTG)

options(scipen=999)

# Exibir os Odds Ratios calculados
print(odds_table_all)

# Criar dataframe dos odds ratios para cada modelo
odds_table_fSFG_full$Model <- "fSFG"
odds_table_fLTG_full$Model <- "fLTG"

# Unir os dois dataframes
odds_table_all <- rbind(odds_table_fSFG_full, odds_table_fLTG_full)

# Criar o gráfico
ggplot(odds_table_all, aes(x = reorder(Variable, OR), y = OR, ymin = Lower95, ymax = Upper95, color = Model)) +
  geom_pointrange(size = 1) +  # Ajustando espessura dos pontos e linhas
  coord_flip() +
  theme_minimal() +  # Mantendo um tema limpo
  labs(x = "", y = "Odds Ratio") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  scale_color_manual(values = c("#e63946", "#4361ee")) +  # Cores personalizadas
  theme_Publication() +  # Substituindo theme_Publication() por um tema similar
  theme(
    text = element_text(family = "sans"),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.position = "bottom",
    legend.spacing.y = unit(1.0, 'cm'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  guides(
    color = guide_legend(nrow = 1, byrow = FALSE)
  )

# Figuras dos modelos ------------------------------------------------------------------------------------------------

if(zmax == 0.03){
  velDisp1   <- 50
  velDisp2   <- 75
}

if(zmax == 0.1){
  velDisp1   <- 115
  velDisp2   <- 150
}

vetor_prob <- c(min(data.s$logvelDisp_e), log10(velDisp1), log10(velDisp2), max(data.s$logvelDisp_e))
#vetor_prob <- c(0.00, 1/3, 2/3, 1.00)

a=1
b=1
c=1

modelos <- data.frame(logvelDisp_e = numeric(),
                      logRproj_rvir = numeric(),
                      
                      pred_fSFG_logito = numeric(),
                      pred_fSFG_logito_se = numeric(),
                      pred_fSFG_prob = numeric(),
                      prob_fSFG_upr = numeric(),
                      prob_fSFG_lwr = numeric(),
                      
                      pred_fLTG_logito = numeric(),
                      pred_fLTG_logito_se = numeric(),
                      pred_fLTG_prob = numeric(),
                      prob_fLTG_upr = numeric(),
                      prob_fLTG_lwr = numeric(),
                      
                      painel_label_logvelDisp_e = character())

# Calculando a fração
bins_total <- data.frame(painel_logvelDisp_e = character(),
                         
                         ngal_painel = numeric(),
                         bin_logRproj_rvir = character(),
                         meio_bin = numeric(),
                         ngal_bin = numeric(),
                         
                         ngal_SF = numeric(),
                         ngal_Q = numeric(),
                         
                         ngal_LT = numeric(),
                         ngal_ET = numeric(),
                         
                         fSFG = numeric(),
                         fSFG_up = numeric(),
                         fSFG_lw = numeric(),
                         
                         fLTG = numeric(),
                         fLTG_up = numeric(),
                         fLTG_lw = numeric())

for (a in 1:(length(vetor_prob)-1)) {
  
  # Probabilidades para quantis:
  prob_lw_logvelDisp_e <- vetor_prob[a]
  prob_up_logvelDisp_e <- vetor_prob[a+1]
  
  # Quantis:
  #quantil_lw_logvelDisp_e <- unname(quantile(data.s$logvelDisp_e, probs = prob_lw_logvelDisp_e))
  #quantil_up_logvelDisp_e <- unname(quantile(data.s$logvelDisp_e, probs = prob_up_logvelDisp_e))
  
  quantil_lw_logvelDisp_e <- prob_lw_logvelDisp_e
  quantil_up_logvelDisp_e <- prob_up_logvelDisp_e
  
  # Definindo galáxias do painel: 
  painel <- 
    data.s$logvelDisp_e >= quantil_lw_logvelDisp_e & 
    data.s$logvelDisp_e < quantil_up_logvelDisp_e 
  
  data.s_painel   <- data.s[painel,]
  ngal_painel <- nrow(data.s_painel)
  
  painel_label_logvelDisp_e <- paste0(
    round(10^quantil_lw_logvelDisp_e, digits = 0), " até ",
    round(10^quantil_up_logvelDisp_e, digits = 0))
  
  painel_label <- paste0("logvelDisp_e = ", painel_label_logvelDisp_e)
  
  # Modelo:
  
  # Densidades: 
  densidade_logvelDisp_e <- density(data.s_painel$logvelDisp_e)
  
  # Modas:
  moda_logvelDisp_e <- densidade_logvelDisp_e$x[which.max(densidade_logvelDisp_e$y)]
  
  # Medianas:
  median_logvelDisp_e <- median(data.s_painel$logvelDisp_e)
  
  # Médias:
  mean_logvelDisp_e <- mean(data.s_painel$logvelDisp_e)
  
  # ------------------------------------
  
  #bins_logRproj_rvir <- unname(quantile(data.s$logRproj_rvir, probs = seq(0, 1, by = 0.05)))
  bins_logRproj_rvir <- seq(min(data.s$logRproj_rvir), max(data.s$logRproj_rvir), by = 0.01)
  
  modelo      <- data.frame(logvelDisp_e = median_logvelDisp_e, logRproj_rvir = bins_logRproj_rvir)
  
  pred_fSFG_logito <- predict(fit_fSFG, type = "link", se.fit = TRUE, newdata = modelo)
  pred_fLTG_logito <- predict(fit_fLTG, type = "link", se.fit = TRUE, newdata = modelo)
  
  modelo$pred_fSFG_logito    <- pred_fSFG_logito$fit
  modelo$pred_fSFG_logito_se <- pred_fSFG_logito$se.fit
  modelo$pred_fSFG_prob      <- fit_fSFG$family$linkinv(pred_fSFG_logito$fit)
  modelo$prob_fSFG_upr       <- fit_fSFG$family$linkinv(pred_fSFG_logito$fit + (critval * pred_fSFG_logito$se.fit))
  modelo$prob_fSFG_lwr       <- fit_fSFG$family$linkinv(pred_fSFG_logito$fit - (critval * pred_fSFG_logito$se.fit))
  
  modelo$pred_fLTG_logito    <- pred_fLTG_logito$fit
  modelo$pred_fLTG_logito_se <- pred_fLTG_logito$se.fit
  modelo$pred_fLTG_prob      <- fit_fLTG$family$linkinv(pred_fLTG_logito$fit)
  modelo$prob_fLTG_upr       <- fit_fLTG$family$linkinv(pred_fLTG_logito$fit + (critval * pred_fLTG_logito$se.fit))
  modelo$prob_fLTG_lwr       <- fit_fLTG$family$linkinv(pred_fLTG_logito$fit - (critval * pred_fLTG_logito$se.fit))
  
  modelo$painel_label_logvelDisp_e <- painel_label_logvelDisp_e
  
  modelos <- rbind(modelos, modelo)
  
  aux <- quantile(data.s_painel$logRproj_rvir, probs = seq(0,1, by = 0.2))
  #aux <- seq(min(data.s$logRproj_rvir), max(data.s$logRproj_rvir), by = 0.25)
  
  for (i in 1:(length(aux)-1)) {
    
    bin <- painel &
      data.s$logRproj_rvir >= aux[i] & 
      data.s$logRproj_rvir < aux[i + 1]
    
    data.s_bin   <- data.s[bin,]
    ngal_bin <- nrow(data.s_bin)
    
    meio_bin <- aux[i] + abs(aux[i+1] - aux[i])/2
    
    aux_label <- paste0(
      round(aux[i], digits = 3), " até ",
      round(aux[i + 1], digits = 3))
    
    ngal_SF <- sum(data.s_bin$SF == 1)
    ngal_Q  <- sum(data.s_bin$SF == 0)
    
    ngal_LT <- sum(data.s_bin$LT == 1)
    ngal_ET <- sum(data.s_bin$LT == 0)
    
    fSFG    <- ngal_SF/ngal_bin
    fSFG_up <- binom.confint(x = ngal_SF, n = ngal_bin, conf.level = 0.95, method = "bayes")$upper
    fSFG_lw <- binom.confint(x = ngal_SF, n = ngal_bin, conf.level = 0.95, method = "bayes")$lower
    
    fLTG    <- ngal_LT/ngal_bin
    fLTG_up <- binom.confint(x = ngal_LT, n = ngal_bin, conf.level = 0.95, method = "bayes")$upper
    fLTG_lw <- binom.confint(x = ngal_LT, n = ngal_bin, conf.level = 0.95, method = "bayes")$lower
    
    bins <- data.frame(painel_logvelDisp_e = painel_label_logvelDisp_e,
                       
                       ngal_painel = ngal_painel,
                       
                       bin_logRproj_rvir = aux_label,
                       meio_bin = meio_bin,
                       ngal_bin = ngal_bin,
                       
                       ngal_SF = ngal_SF,
                       ngal_Q = ngal_Q,
                       
                       ngal_LT = ngal_LT,
                       ngal_ET = ngal_ET,
                       
                       fSFG = fSFG,
                       fSFG_up = fSFG_up,
                       fSFG_lw = fSFG_lw,
                       
                       fLTG = fLTG,
                       fLTG_up = fLTG_up,
                       fLTG_lw = fLTG_lw)
    
    bins_total <- rbind(bins_total, bins)
  }
}

colnames(modelos)[13] <- "painel_logvelDisp_e"

lv1 <- paste0(round(min(data.s$velDisp_e), digits = 0), " até ", velDisp1)
lv2 <- paste0(velDisp1, " até ", velDisp2)
lv3 <- paste0(velDisp2, " até ", round(max(data.s$velDisp_e), digits = 0))

levels <- c(lv1, lv2, lv3)

modelos$painel_logvelDisp_e    <- factor(modelos$painel_logvelDisp_e, levels = levels)
bins_total$painel_logvelDisp_e <- factor(bins_total$painel_logvelDisp_e, levels = levels)

colors <- c("fSFG" = "blue", "fLTG" = "red")

ggplot() + 
  geom_line(data = modelos, aes(x = logRproj_rvir, y = pred_fSFG_prob), linewidth = 1, linetype = "dotdash", color = '#4361ee') + 
  geom_ribbon(data = modelos, aes(ymin = prob_fSFG_lwr, ymax = prob_fSFG_upr, x = logRproj_rvir), alpha = 0.5, inherit.aes = FALSE, fill = 'blue') + 
  
  geom_line(data = modelos, aes(x = logRproj_rvir, y = pred_fLTG_prob), linewidth = 1, linetype = "dotdash", color = '#e63946') + 
  geom_ribbon(data = modelos, aes(ymin = prob_fLTG_lwr, ymax = prob_fLTG_upr, x = logRproj_rvir), alpha = 0.5, inherit.aes = FALSE, fill = 'red') + 
  
  geom_point(data = bins_total, aes(x = meio_bin, y = fSFG), size = 2, color = '#4361ee') +
  geom_errorbar(data = bins_total, aes(ymin = fSFG_lw, ymax = fSFG_up, x = meio_bin), width = 0.1, color = '#4361ee') + 
  
  geom_point(data = bins_total, aes(x = meio_bin, y = fLTG), size = 2, color = '#e63946') +
  geom_errorbar(data = bins_total, aes(ymin = fLTG_lw, ymax = fLTG_up, x = meio_bin), width = 0.1, color = '#e63946') + 
  
  facet_grid(. ~ painel_logvelDisp_e) + 
  
  ylab("Fraction") + 
  
  scale_x_continuous(sec.axis = sec_axis(~ . , name = label_logvelDisp_e, breaks = NULL, labels = NULL)) + 
  
  scale_color_manual(values = c("#e63946", "#4361ee")) + 
  scale_fill_manual(values = c("#e63946", "#4361ee")) + 
  
  theme_Publication() + 
  theme(text = element_text(family = "sans"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.title.y.right = element_text(size = 14),
        axis.title.x.top = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.position = "bottom",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) + 
  
  guides(fill = guide_legend(nrow = 1, byrow = F),
         color = guide_legend(nrow = 1, byrow = F),
         legend.spacing.y = unit(1.0, 'cm'))


