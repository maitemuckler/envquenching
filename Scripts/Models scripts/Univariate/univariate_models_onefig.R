library(data.table)
library(InformationValue)
library(ggplot2)
library(scales)
library(dplyr)
library(binom)
library(patchwork)
library(purrr)

## Lendo meus códigos ----
source("~/Work/Research/Astronomy/Projects/environmental-quenching/Scripts/Themes/my_theme.R")
source("~/Work/Research/Astronomy/Projects/environmental-quenching/Scripts/Themes/ggplot_theme_Publication-2.R")

figswd <- "~/Work/Research/Astronomy/Projects/environmental-quenching/Figures/Models/univariate_models"
critval <- 1.96 

# Lista de valores para 'preditora'
preditora_list <- c("TType", "logvelDisp_e", "logMstar", "logRproj_rvir")

# Função para gerar gráficos baseados em cada preditora
generate_plot <- function(preditora) {
  
  # Configurações para a variável 'preditora'
  if(preditora == "logvelDisp_e"){preditora_label <- label_logvelDisp_e; n_bins <- 10; binagem <- "cut"; xlim_label <- -Inf}
  if(preditora == "logMstar"){preditora_label <- label_logMstar; n_bins <- 10; binagem <- "cut"; xlim_label <- -Inf}
  if(preditora == "logRproj_rvir"){preditora_label <- label_logRproj_rvir; n_bins <- 10; binagem <- "cut"; xlim_label <- Inf}
  if(preditora == "TType"){preditora_label <- label_TType; n_bins <- 10; binagem <- "cut"; xlim_label <- Inf}
  
  input_data003 <- paste0("inputdata_zmax0.03_Rlim2.5_Ma12.3_flag_good==1_MANGLE_logMstar_min10.5.csv")
  input_data01  <- paste0("inputdata_zmax0.1_Rlim2.5_Ma12.3_flag_good==1_MANGLE_logMstar_min10.5.csv")
  
  colunas <- c("SF_GSWLC", "type", preditora)
  
  df003 <- fread(paste0("~/Work/Research/Astronomy/Data/EnvQuenching/inputModel/GSWLC/", input_data003), select = colunas)
  df01  <- fread(paste0("~/Work/Research/Astronomy/Data/EnvQuenching/inputModel/GSWLC/", input_data01), select = colunas)
  
  colnames(df003)[which(colnames(df003) == preditora)] <- "preditora"
  colnames(df01)[which(colnames(df01) == preditora)] <- "preditora"
  
  # Filtrando apenas Satélites e convertendo variáveis
  df003 <- subset(df003, df003$type == "Satellite")
  df01  <- subset(df01, df01$type == "Satellite")
  
  df003$zmax <- 0.03
  df01$zmax  <- 0.1
  
  df003$SF <- as.factor(ifelse(df003$SF_GSWLC == "Star-forming", 1, 0))
  df01$SF <- as.factor(ifelse(df01$SF_GSWLC == "Star-forming", 1, 0))
  
  # Modelos e previsões
  modelo003 <- glm(SF ~ preditora, family = binomial(link = "logit"), data = df003)
  modelo01  <- glm(SF ~ preditora, family = binomial(link = "logit"), data = df01)
  df003$pred <- predict(modelo003, type = "response")
  df01$pred  <- predict(modelo01, type = "response")
  
  optimal003 <- optimalCutoff(df003$SF, df003$pred)[1]
  optimal01  <- optimalCutoff(df01$SF, df01$pred)[1]
  
  classe_predita003 <- ifelse(df003$pred > optimal003, 1, 0)
  classe_predita01  <- ifelse(df01$pred > optimal01, 1, 0)
  
  # Calcular a taxa de acerto
  acertos003     <- sum(classe_predita003 == df003$SF)
  taxa_acerto003 <- acertos003 / nrow(df003)
  taxa_acerto003 <- round(taxa_acerto003*100, digits = 1)
  
  acertos01     <- sum(classe_predita01 == df01$SF)
  taxa_acerto01 <- acertos01 / nrow(df01)
  taxa_acerto01 <- round(taxa_acerto01*100, digits = 1)
  
  # Previsão para a faixa de valores
  preddata003 <- data.frame(preditora = seq(min(df003$preditora), max(df003$preditora), length = 100))
  preddata01  <- data.frame(preditora = seq(min(df01$preditora), max(df01$preditora), length = 100))
  preds003 <- predict(modelo003, newdata = preddata003, type = "link", se.fit = TRUE)
  preds01  <- predict(modelo01, newdata = preddata01, type = "link", se.fit = TRUE)
  
  preddata003$accuracy <- taxa_acerto003
  preddata01$accuracy  <- taxa_acerto01
  
  # Confiança
  preddata003$lwr <- modelo003$family$linkinv(preds003$fit - critval * preds003$se.fit)
  preddata003$upr <- modelo003$family$linkinv(preds003$fit + critval * preds003$se.fit)
  preddata01$lwr <- modelo01$family$linkinv(preds01$fit - critval * preds01$se.fit)
  preddata01$upr <- modelo01$family$linkinv(preds01$fit + critval * preds01$se.fit)
  
  # Adicionando a coluna zmax
  preddata003$zmax <- 0.03
  preddata01$zmax <- 0.1
  
  # Preparação dos dados de binagem
  bins_temp003 <- df003 %>%
    mutate(bins = if(binagem == "cut") cut(preditora, n_bins) else ntile(preditora, n_bins)) %>%
    group_by(bins) %>%
    summarize(
      mediana_bins = median(preditora),
      proporcao_SF_igual_1 = mean(SF == 1),
      confint_up = binom.confint(sum(SF == 1), n(), conf.level = 0.95, method = "bayes")$upper,
      confint_lw = binom.confint(sum(SF == 1), n(), conf.level = 0.95, method = "bayes")$lower,
      .groups = 'drop'
    )
  bins_temp01 <- df01 %>%
    mutate(bins = if(binagem == "cut") cut(preditora, n_bins) else ntile(preditora, n_bins)) %>%
    group_by(bins) %>%
    summarize(
      mediana_bins = median(preditora),
      proporcao_SF_igual_1 = mean(SF == 1),
      confint_up = binom.confint(sum(SF == 1), n(), conf.level = 0.95, method = "bayes")$upper,
      confint_lw = binom.confint(sum(SF == 1), n(), conf.level = 0.95, method = "bayes")$lower,
      .groups = 'drop'
    )
  
  # Adicionando a coluna zmax antes de combinar
  bins_temp003$zmax <- 0.03
  bins_temp01$zmax <- 0.1
  
  # Combinando dados
  df <- bind_rows(df003, df01)
  bins_temp <- bind_rows(bins_temp003, bins_temp01)
  preddata <- bind_rows(preddata003, preddata01)
  
  # Transformando zmax em fator para usá-lo nos gráficos
  df$zmax <- as.factor(df$zmax)
  bins_temp$zmax <- as.factor(bins_temp$zmax)
  preddata$zmax <- as.factor(preddata$zmax)
  
  # Construindo o gráfico
  g <- ggplot() + 
    geom_line(data = df, aes(x = preditora, y = pred, color = zmax)) +         
    geom_point(data = bins_temp, aes(y = proporcao_SF_igual_1, x = mediana_bins, color = zmax), size = 2) +
    geom_errorbar(data = bins_temp, aes(ymin = confint_lw, ymax = confint_up, x = mediana_bins, color = zmax), width = 0.025, alpha = 0.4) + 
    geom_ribbon(data = preddata, aes(ymin = lwr, ymax = upr, x = preditora, fill = zmax), alpha = 0.3) +
    
    annotate("text", x = -Inf, y = 0, label = "Accuracy: ", hjust = -0.2, vjust = -1, family = "sans", size = 5) + 
    
    geom_label(aes(x = -Inf, y = 0, label = paste(round(taxa_acerto003, 1), "%")), 
               data = data.frame(zmax = factor(0.03)), hjust = -1.5, vjust = 0.1, color = "#8ac926", fill = "white", size = 5) +
    geom_label(aes(x = -Inf, y = 0, label = paste(round(taxa_acerto01, 1), "%")), 
               data = data.frame(zmax = factor(0.1)), hjust = -1.5, vjust = -1, color = "#ffb703", fill = "white", size = 5) +
    
    scale_color_manual(values = c("#8ac926", "#ffb703")) + 
    scale_fill_manual(values = c("#8ac926", "#ffb703")) + 
    
    scale_y_continuous(breaks = pretty_breaks(n = 5)) +
    scale_x_continuous(breaks = pretty_breaks(n = 5)) +
    
    labs(x = preditora_label, y = "fSFG", color = expression(z[max]), fill = expression(z[max])) + 
    coord_cartesian(ylim = c(0, 1)) + 
    
    theme_Publication() + 
    theme(axis.text = element_text(size = 18),
          axis.title = element_text(size = 22),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 22))
  
  return(g)
}

# Gerando o gráfico para cada preditora e combinando com patchwork
plots <- map(preditora_list, generate_plot)
combined_plot <- wrap_plots(plots) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
combined_plot

#Salvando a figura combinada
ggsave(path = figswd, plot = combined_plot,
       filename = paste0("univariate_models.pdf"),
       device = cairo_pdf, width = width_figs, height = height_figs, units = "in", dpi = 600)
