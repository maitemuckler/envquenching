# Diretórios ----
setwd("~/Work/Research/Astronomy/Projects/envquenching/")
wddata <- "~/Work/Research/Astronomy/Data/"
figswd <- "~/Work/Research/Astronomy/Projects/envquenching/Figures/Letter/"

# Bibliotecas ----
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

# Códigos extras ----
source("Scripts/Themes/my_theme.R")

# Opções ----
options(scipen = 999)

critval   <- 1.96 
TType_lim <- 0

width_figs  <- 9
height_figs <- 5

# Dados ----
input_data <- paste0("/inputdata_zmax0.1_Rlim2.5_Ma12.3_flag_good==1_MANGLE_logMstar_min10.5.csv")
df         <- fread(paste0(wddata, "inputModel/GSWLC", input_data))

data.s <- subset(df, df$type == "Satellite")
data.c <- subset(df, df$type == "Central")

# Modelo para fSFG ----

data.s$SF  <- ifelse(data.s$SF_GSWLC == "Star-forming", 1, 0)
fit_fSFG  <- glm(SF ~ logvelDisp_e + logRproj_rvir, family = binomial(link = "logit"), data = data.s)

# Resumo do modelo
summary(fit_fSFG)

# Intervalo de confiança para os coeficientes
confint(fit_fSFG)

# Odds Ratio
OR_fSFG <- exp(coef(fit_fSFG))
OR_fSFG

# IC para Odds Ratio
conf_int_OR_fSFG <- exp(confint(fit_fSFG))
conf_int_OR_fSFG

# Tabela para OR
OR_table_fSFG <- data.frame(
  Variable = names(OR_fSFG),
  OR = OR_fSFG,
  Lower95 = conf_int_OR_fSFG[,1],
  Upper95 = conf_int_OR_fSFG[,2]
)

rownames(OR_table_fSFG) <- NULL

OR_table_fSFG

# Importância das covariáveis
importance_fSFG <- as.data.frame(varImp(fit_fSFG))
importance_fSFG <- data.frame(overall = importance_fSFG$Overall, names = rownames(importance_fSFG))
importance_fSFG[order(importance_fSFG$overall, decreasing = T),]

# Predição
data.s$pred_fSFG <- predict(fit_fSFG, data.s, type = "response")
optimal_fSFG     <- optimalCutoff(data.s$SF, data.s$pred_fSFG)[1]
misClass_fSFG    <- misClassError(data.s$SF, data.s$pred_fSFG, threshold = optimal_fSFG)
acuracia_fSFG    <- 1 - misClass_fSFG
acuracia_fSFG

# Matriz de confusão
data.s$pred_class_fSFG <- ifelse(data.s$pred_fSFG >= optimal_fSFG, "Star-forming", "Quiescent") # Converter probabilidades em classes
conf_matrix_fSFG       <- confusionMatrix(factor(data.s$pred_class_fSFG), factor(data.s$SF_GSWLC))
conf_matrix_fSFG

# Extraindo métricas de desempenho para fSFG
acc_fSFG         <- conf_matrix_fSFG$overall["Accuracy"]
sensitivity_fSFG <- conf_matrix_fSFG$byClass["Sensitivity"]
specificity_fSFG <- conf_matrix_fSFG$byClass["Specificity"]

# Curva ROC para fSFG
roc_fSFG <- roc(data.s$SF, data.s$pred_fSFG)
plot(roc_fSFG, main="ROC Curve - fSFG", col="blue", lwd=2)
cat("AUC para fSFG:", auc(roc_fSFG), "\n")

# Modelo para fLTG ----

data.s$LT <- ifelse(data.s$TType >= TType_lim, 1, 0)
data.s$LT <- as.factor(data.s$LT)

fit_fLTG <- glm(LT ~ logvelDisp_e + logRproj_rvir, family = binomial(link = "logit"), data = data.s)

# Resumo do modelo
summary(fit_fLTG)

# Intervalo de confiança para os coeficientes
confint(fit_fLTG)

# Odds Ratio
OR_fLTG <- exp(coef(fit_fLTG))
OR_fLTG

# IC para Odds Ratio
conf_int_OR_fLTG <- exp(confint(fit_fLTG))
conf_int_OR_fLTG

# Tabela para OR
OR_table_fLTG <- data.frame(
  Variable = names(OR_fLTG),
  OR = OR_fLTG,
  Lower95 = conf_int_OR_fLTG[,1],
  Upper95 = conf_int_OR_fLTG[,2]
)

rownames(OR_table_fLTG) <- NULL

OR_table_fLTG

# Importância das covariáveis
importance_fLTG <- as.data.frame(varImp(fit_fLTG))
importance_fLTG <- data.frame(overall = importance_fLTG$Overall, names = rownames(importance_fLTG))
importance_fLTG[order(importance_fLTG$overall, decreasing = T),]

# Predição
data.s$pred_fLTG <- predict(fit_fLTG, data.s, type = "response")
optimal_fLTG     <- optimalCutoff(data.s$LT, data.s$pred_fLTG)[1]
misClass_fLTG    <- misClassError(data.s$LT, data.s$pred_fLTG, threshold = optimal_fLTG)
acuracia_fLTG    <- 1 - misClass_fLTG
acuracia_fLTG

# Matriz de confusão
data.s$pred_class_fLTG <- ifelse(data.s$pred_fLTG >= optimal_fLTG, "Late-type", "Early-type") # Converter probabilidades em classes
conf_matrix_fLTG       <- confusionMatrix(factor(data.s$pred_class_fLTG), factor(data.s$morph_char))
conf_matrix_fLTG

# Extraindo métricas de desempenho para fSFG
acc_fLTG         <- conf_matrix_fLTG$overall["Accuracy"]
sensitivity_fLTG <- conf_matrix_fLTG$byClass["Sensitivity"]
specificity_fLTG <- conf_matrix_fLTG$byClass["Specificity"]

# Curva ROC para fSFG
roc_fLTG <- roc(data.s$LT, data.s$pred_fLTG)
plot(roc_fLTG, main="ROC Curve - fLTG", col="blue", lwd=2)
cat("AUC para fLTG:", auc(roc_fLTG), "\n")

# --------------------------------------------------------------------------------------------------------

# Código para Criar o Gráfico de Odds Ratios

# Adicionar rótulos para identificar os modelos
OR_table_fSFG$Model <- "fSFG"
OR_table_fLTG$Model <- "fLTG"

# Unir os dois dataframes
OR_table <- rbind(OR_table_fSFG[-1,], OR_table_fLTG[-1,])
OR_table <- rbind(OR_table_fSFG, OR_table_fLTG)
rownames(OR_table) <- NULL
OR_table

# Criar o gráfico
ggplot(OR_table, aes(x = reorder(Variable, OR), y = OR, ymin = Lower95, ymax = Upper95, color = Model)) +
  geom_pointrange(size = 1) +  # Ajustando espessura dos pontos e linhas
  coord_flip() +
  facet_grid(. ~ Variable, scales = "free_x") + 
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
# --------------------------------------------------------------------------------------------------------

## FIGURAS LOGÍSTICAS:

# Definição dos quantis
velDisp1 <- 120
velDisp2 <- 153
vetor_prob <- c(min(data.s$logvelDisp_e), log10(velDisp1), log10(velDisp2), max(data.s$logvelDisp_e))

#vetor_prob <- quantile(data.s$logvelDisp_e, probs = seq(0, 1, length.out = 4))

# Listas para otimizar a construção dos data frames
modelos_list    <- list()
bins_total_list <- list()

# Loop para calcular os modelos e bins
for (a in 1:(length(vetor_prob)-1)) {
  
  quantil_lw_logvelDisp_e <- vetor_prob[a]
  quantil_up_logvelDisp_e <- vetor_prob[a+1]
  
  painel <- data.s$logvelDisp_e >= quantil_lw_logvelDisp_e & 
    data.s$logvelDisp_e < quantil_up_logvelDisp_e 
  
  data.s_painel <- data.s[painel,]
  ngal_painel   <- nrow(data.s_painel)
  
  painel_logvelDisp_e <- paste0("[",
                                round(10^quantil_lw_logvelDisp_e, digits = 0), ", ",
                                round(10^quantil_up_logvelDisp_e, digits = 0), ")")
  
  # Criando os modelos
  bins_logRproj_rvir <- seq(min(data.s_painel$logRproj_rvir), max(data.s_painel$logRproj_rvir), by = 0.01)
  
  modelo <- data.frame(logvelDisp_e = median(data.s_painel$logvelDisp_e), 
                       logRproj_rvir = bins_logRproj_rvir)
  
  pred_fSFG <- predict(fit_fSFG, type = "link", se.fit = TRUE, newdata = modelo)
  pred_fLTG <- predict(fit_fLTG, type = "link", se.fit = TRUE, newdata = modelo)
  
  modelo$pred_fSFG_prob <- fit_fSFG$family$linkinv(pred_fSFG$fit)
  modelo$prob_fSFG_upr  <- fit_fSFG$family$linkinv(pred_fSFG$fit + (critval * pred_fSFG$se.fit))
  modelo$prob_fSFG_lwr  <- fit_fSFG$family$linkinv(pred_fSFG$fit - (critval * pred_fSFG$se.fit))
  
  modelo$pred_fLTG_prob <- fit_fLTG$family$linkinv(pred_fLTG$fit)
  modelo$prob_fLTG_upr  <- fit_fLTG$family$linkinv(pred_fLTG$fit + (critval * pred_fLTG$se.fit))
  modelo$prob_fLTG_lwr  <- fit_fLTG$family$linkinv(pred_fLTG$fit - (critval * pred_fLTG$se.fit))
  
  modelo$painel_logvelDisp_e <- painel_logvelDisp_e
  modelos_list[[a]] <- modelo
  
  # Criando os bins
  aux <- quantile(data.s_painel$logRproj_rvir, probs = seq(0, 1, by = 0.2))
  
  for (i in 1:(length(aux)-1)) {
    
    bin <- painel & data.s$logRproj_rvir >= aux[i] & data.s$logRproj_rvir < aux[i+1]
    data.s_bin <- data.s[bin,]
    ngal_bin <- nrow(data.s_bin)
    
    if (ngal_bin > 0) {
      meio_bin <- mean(c(aux[i], aux[i+1]))
      
      ngal_SF <- sum(data.s_bin$SF == 1)
      ngal_LT <- sum(data.s_bin$LT == 1)
      
      fSFG <- ngal_SF / ngal_bin
      fLTG <- ngal_LT / ngal_bin
      
      fSFG_up <- binom.confint(ngal_SF, ngal_bin, conf.level = 0.95, method = "bayes")$upper
      fSFG_lw <- binom.confint(ngal_SF, ngal_bin, conf.level = 0.95, method = "bayes")$lower
      
      fLTG_up <- binom.confint(ngal_LT, ngal_bin, conf.level = 0.95, method = "bayes")$upper
      fLTG_lw <- binom.confint(ngal_LT, ngal_bin, conf.level = 0.95, method = "bayes")$lower
      
      bins_total_list[[length(bins_total_list) + 1]] <- data.frame(
        painel_logvelDisp_e = painel_logvelDisp_e,
        meio_bin = meio_bin,
        fSFG = fSFG, fSFG_up = fSFG_up, fSFG_lw = fSFG_lw,
        fLTG = fLTG, fLTG_up = fLTG_up, fLTG_lw = fLTG_lw)
    }
  }
}

# Probs específicas: Radial variation of the fraction ~ delta fraction

# Criando um data frame vazio para armazenar os resultados
resultados_df <- data.frame(
  logvelDisp_e = numeric(),
  pred_fSFG_prob_min = numeric(),
  pred_fLTG_prob_min = numeric(),
  pred_fSFG_prob_max = numeric(),
  pred_fLTG_prob_max = numeric(),
  pred_fSFG_prob_abs = numeric(),
  pred_fLTG_prob_abs = numeric(),
  stringsAsFactors = FALSE
)

# Obtendo os valores únicos de logvelDisp_e
valores_logvelDisp_e <- unique(modelos$logvelDisp_e)

# Loop sobre os três painéis
for (i in seq_along(valores_logvelDisp_e)) {
  painel <- subset(modelos, modelos$logvelDisp_e == valores_logvelDisp_e[i])
  
  # Obtendo os índices dos valores desejados
  min_logRproj_idx <- which.min(painel$logRproj_rvir)
  max_logRproj_idx <- which.max(painel$logRproj_rvir)
  min_abs_logRproj_idx <- which.min(abs(painel$logRproj_rvir))
  
  # Criando uma nova linha de resultados
  nova_linha <- data.frame(
    logvelDisp_e = valores_logvelDisp_e[i],
    
    # Valores extremos
    pred_fSFG_prob_min = painel$pred_fSFG_prob[min_logRproj_idx],
    pred_fLTG_prob_min = painel$pred_fLTG_prob[min_logRproj_idx],
    
    pred_fSFG_prob_max = painel$pred_fSFG_prob[max_logRproj_idx],
    pred_fLTG_prob_max = painel$pred_fLTG_prob[max_logRproj_idx],
    
    pred_fSFG_prob_abs = painel$pred_fSFG_prob[min_abs_logRproj_idx],
    pred_fLTG_prob_abs = painel$pred_fLTG_prob[min_abs_logRproj_idx]
    
  )
  
  # Adicionando a nova linha ao data frame de resultados
  resultados_df <- rbind(resultados_df, nova_linha)
}

# Diferença para fSFG:
round(((resultados_df$pred_fSFG_prob_max - resultados_df$pred_fSFG_prob_min) * 100), digits = 2)

# Diferença para fLTG:
round(((resultados_df$pred_fLTG_prob_max - resultados_df$pred_fLTG_prob_min) * 100), digits = 2)

# Convertendo listas para data frames
modelos    <- do.call(rbind, modelos_list)
bins_total <- do.call(rbind, bins_total_list)

# Garantindo que painel_logvelDisp_e é um fator
modelos$painel_logvelDisp_e <- factor(modelos$painel_logvelDisp_e,
                                      levels = c(paste0("[", round(10^vetor_prob)[1], ", ", round(10^vetor_prob)[2], ")"), 
                                                 paste0("[", round(10^vetor_prob)[2], ", ", round(10^vetor_prob)[3], ")"),
                                                 paste0("[", round(10^vetor_prob)[3], ", ", round(10^vetor_prob)[4], ")")))
bins_total$painel_logvelDisp_e <- factor(bins_total$painel_logvelDisp_e,
                                         levels = c(paste0("[", round(10^vetor_prob)[1], ", ", round(10^vetor_prob)[2], ")"), 
                                                    paste0("[", round(10^vetor_prob)[2], ", ", round(10^vetor_prob)[3], ")"),
                                                    paste0("[", round(10^vetor_prob)[3], ", ", round(10^vetor_prob)[4], ")")))

# Plotagem
ggplot() + 
  geom_line(data = modelos, aes(x = logRproj_rvir, y = pred_fSFG_prob, color = "fSFG", group = painel_logvelDisp_e), linewidth = 1) + 
  geom_ribbon(data = modelos, aes(ymin = prob_fSFG_lwr, ymax = prob_fSFG_upr, x = logRproj_rvir, fill = "fSFG", group = painel_logvelDisp_e), alpha = 0.5) + 
  
  geom_line(data = modelos, aes(x = logRproj_rvir, y = pred_fLTG_prob, color = "fLTG", group = painel_logvelDisp_e), linewidth = 1) + 
  geom_ribbon(data = modelos, aes(ymin = prob_fLTG_lwr, ymax = prob_fLTG_upr, x = logRproj_rvir, fill = "fLTG", group = painel_logvelDisp_e), alpha = 0.5) + 
  
  geom_point(data = bins_total, aes(x = meio_bin, y = fSFG, color = "fSFG", group = painel_logvelDisp_e), size = 2) +
  geom_errorbar(data = bins_total, aes(ymin = fSFG_lw, ymax = fSFG_up, x = meio_bin, color = "fSFG", group = painel_logvelDisp_e), width = 0.1) + 
  
  geom_point(data = bins_total, aes(x = meio_bin, y = fLTG, color = "fLTG", group = painel_logvelDisp_e), size = 2) +
  geom_errorbar(data = bins_total, aes(ymin = fLTG_lw, ymax = fLTG_up, x = meio_bin, color = "fLTG", group = painel_logvelDisp_e), width = 0.1) + 
  
  facet_grid(. ~ painel_logvelDisp_e) + 
  
  ylab("Fraction") + 
  xlab(label_logRproj_rvir) + 
  
  scale_color_manual(values = c("fSFG" = "#9B5DE5", "fLTG" = "#7FD23F"), 
                     name = "",
                     labels = c("fSFG" = expression(f[SFG]), "fLTG" = expression(f[LTG])),
                     guide = guide_legend(override.aes = list(linewidth = 1.5, size = 4))) + 
  scale_fill_manual(values = c("fSFG" = "#9B5DE5", "fLTG" = "#7FD23F")) + 
  
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  scale_x_continuous(sec.axis = sec_axis(~ . , expression(paste(sigma[e], " (km/s)")), breaks = NULL, labels = NULL)) + 
  
  theme_Publication() + 
  theme(text = element_text(size = 18),
        legend.text = element_text(size = 20),
        legend.position = "inside",
        legend.position.inside = c(0.16, 0.11),
        legend.key.size = unit(1, "cm")) + 
  guides(fill = "none")

ggsave(path = figswd,
       filename = paste0("logistic_z010.pdf"),
       device = cairo_pdf, width = width_figs, height = height_figs, 
       units = "in", dpi = 600)
