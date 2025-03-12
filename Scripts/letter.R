setwd("~/Work/Research/Astronomy/Projects/envquenching/")

library(data.table)
library(binom)
library(ggplot2)
library(InformationValue)
library(caret)
library(dplyr)
library(scales)
library(viridis)
library(ggthemes)
library(pROC)

options(scipen = 999)

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

# Modelos ---------------------------------------------------------------------------------------

# fSFG zmax = 0.03

fit_fSFG <- glm(SF ~ logvelDisp_e + logRproj_rvir, family = binomial(link = "logit"), data = data.s)

summary(fit_fSFG)
confint(fit_fSFG)

OR_fSFG <- exp(coef(fit_fSFG))

importance_fSFG <- as.data.frame(varImp(fit_fSFG))
importance_fSFG <- data.frame(overall = importance_fSFG$Overall, names = rownames(importance_fSFG))
importance_fSFG[order(importance_fSFG$overall, decreasing = T),]

# Predição
data.s$pred_fSFG <- predict(fit_fSFG, data.s, type = "response")
optimal_fSFG     <- optimalCutoff(data.s$SF, data.s$pred_fSFG)[1]
misClass_fSFG    <- misClassError(data.s$SF, data.s$pred_fSFG, threshold = optimal_fSFG)
acuracia_fSFG    <- 1 - misClass_fSFG

# Matriz de confusão
 
data.s$pred_class_fSFG <- ifelse(data.s$pred_fSFG >= optimal_fSFG, "Star-forming", "Quiescent") # Converter probabilidades em classes
conf_matrix_fSFG       <- confusionMatrix(factor(data.s$pred_class_fSFG), factor(data.s$SF_char))
conf_matrix_fSFG

# Extraindo métricas de desempenho para fSFG
acc_fSFG         <- conf_matrix_fSFG$overall["Accuracy"]
sensitivity_fSFG <- conf_matrix_fSFG$byClass["Sensitivity"]
specificity_fSFG <- conf_matrix_fSFG$byClass["Specificity"]

# Curva ROC para fSFG
roc_fSFG <- roc(data.s$SF, data.s$pred_fSFG)
plot(roc_fSFG, main="ROC Curve - fSFG", col="blue", lwd=2)
cat("AUC para fSFG:", auc(roc_fSFG), "\n")

# fLTG zmax = 0.03

fit_fLTG <- glm(LT ~ logvelDisp_e + logRproj_rvir, family = binomial(link = "logit"), data = data.s)

summary(fit_fLTG)
confint(fit_fLTG)

OR_fSFG <- exp(coef(fit_fLTG))

importance_fLTG <- as.data.frame(varImp(fit_fLTG))
importance_fLTG <- data.frame(overall = importance_fLTG$Overall, names = rownames(importance_fLTG))
importance_fLTG[order(importance_fLTG$overall, decreasing = T),]

# Predição
data.s$pred_fLTG <- predict(fit_fLTG, data.s, type = "response")
optimal_fLTG     <- optimalCutoff(data.s$LT, data.s$pred_fLTG)[1]
misClass_fLTG    <- misClassError(data.s$LT, data.s$pred_fLTG, threshold = optimal_fLTG)
acuracia_fLTG    <- 1 - misClass_fLTG

# Matriz de confusão

data.s$pred_class_fLTG <- ifelse(data.s$pred_fLTG >= optimal_fLTG, "Late-type", "Early-type") # Converter probabilidades em classes
conf_matrix_fLTG       <- confusionMatrix(factor(data.s$pred_class_fLTG), factor(data.s$morph_char))
conf_matrix_fLTG

# Extraindo métricas de desempenho para fLTG
acc_fLTG         <- conf_matrix_fLTG$overall["Accuracy"]
sensitivity_fLTG <- conf_matrix_fLTG$byClass["Sensitivity"]
specificity_fLTG <- conf_matrix_fLTG$byClass["Specificity"]

# Curva ROC para fLTG
roc_fLTG <- roc(data.s$LT, data.s$pred_fLTG)
plot(roc_fLTG, main="ROC Curve - fLTG", col="red", lwd=2)
cat("AUC para fLTG:", auc(roc_fLTG), "\n")

# Modelos ---------------------------------------------------------------------------------------

# fSFG zmax = 0.03

fit_fSFG <- glm(SF ~ logvelDisp_e + logRproj_rvir, family = binomial(link = "logit"), data = data.s)

summary(fit_fSFG)
confint(fit_fSFG)

OR_fSFG <- exp(coef(fit_fSFG))

importance_fSFG <- as.data.frame(varImp(fit_fSFG))
importance_fSFG <- data.frame(overall = importance_fSFG$Overall, names = rownames(importance_fSFG))
importance_fSFG[order(importance_fSFG$overall, decreasing = T),]

# Predição
data.s$pred_fSFG <- predict(fit_fSFG, data.s, type = "response")
optimal_fSFG     <- optimalCutoff(data.s$SF, data.s$pred_fSFG)[1]
misClass_fSFG    <- misClassError(data.s$SF, data.s$pred_fSFG, threshold = optimal_fSFG)
acuracia_fSFG    <- 1 - misClass_fSFG

# Matriz de confusão
 
data.s$pred_class_fSFG <- ifelse(data.s$pred_fSFG >= optimal_fSFG, "Star-forming", "Quiescent") # Converter probabilidades em classes
conf_matrix_fSFG       <- confusionMatrix(factor(data.s$pred_class_fSFG), factor(data.s$SF_char))
conf_matrix_fSFG

# Extraindo métricas de desempenho para fSFG
acc_fSFG         <- conf_matrix_fSFG$overall["Accuracy"]
sensitivity_fSFG <- conf_matrix_fSFG$byClass["Sensitivity"]
specificity_fSFG <- conf_matrix_fSFG$byClass["Specificity"]

# Curva ROC para fSFG
roc_fSFG <- roc(data.s$SF, data.s$pred_fSFG)
plot(roc_fSFG, main="ROC Curve - fSFG", col="blue", lwd=2)
cat("AUC para fSFG:", auc(roc_fSFG), "\n")

# fLTG zmax = 0.03

fit_fLTG <- glm(LT ~ logvelDisp_e + logRproj_rvir, family = binomial(link = "logit"), data = data.s)

summary(fit_fLTG)
confint(fit_fLTG)

OR_fSFG <- exp(coef(fit_fLTG))

importance_fLTG <- as.data.frame(varImp(fit_fLTG))
importance_fLTG <- data.frame(overall = importance_fLTG$Overall, names = rownames(importance_fLTG))
importance_fLTG[order(importance_fLTG$overall, decreasing = T),]

# Predição
data.s$pred_fLTG <- predict(fit_fLTG, data.s, type = "response")
optimal_fLTG     <- optimalCutoff(data.s$LT, data.s$pred_fLTG)[1]
misClass_fLTG    <- misClassError(data.s$LT, data.s$pred_fLTG, threshold = optimal_fLTG)
acuracia_fLTG    <- 1 - misClass_fLTG

# Matriz de confusão

data.s$pred_class_fLTG <- ifelse(data.s$pred_fLTG >= optimal_fLTG, "Late-type", "Early-type") # Converter probabilidades em classes
conf_matrix_fLTG       <- confusionMatrix(factor(data.s$pred_class_fLTG), factor(data.s$morph_char))
conf_matrix_fLTG

# Extraindo métricas de desempenho para fLTG
acc_fLTG         <- conf_matrix_fLTG$overall["Accuracy"]
sensitivity_fLTG <- conf_matrix_fLTG$byClass["Sensitivity"]
specificity_fLTG <- conf_matrix_fLTG$byClass["Specificity"]

# Curva ROC para fLTG
roc_fLTG <- roc(data.s$LT, data.s$pred_fLTG)
plot(roc_fLTG, main="ROC Curve - fLTG", col="red", lwd=2)
cat("AUC para fLTG:", auc(roc_fLTG), "\n")
