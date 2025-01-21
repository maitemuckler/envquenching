setwd("~/Work/Research/Astronomy/Projects/environmental-quenching/Scripts/Models scripts/")

library(data.table)
library(binom)
library(ggplot2)
library(InformationValue)
library(caret)
library(scales)

source("~/Work/Research/Astronomy/Projects/environmental-quenching/Scripts/Themes/my_theme.R")

zmax       <- 0.1
TType_lim  <- 0
critval    <- 1.96 

wddata     <- "~/Work/Research/Astronomy/Data/environmental-quenching-data/"
input_data <- paste0("inputdata_zmax",zmax,"_Rlim2.5_Ma12.3_flag_good==1_MANGLE_logMstar_min10.5.csv")

df    <- fread(paste0(wddata, "inputModel/GSWLC/", input_data))
df    <- subset(df, df$type == "Satellite")

df$SF <- ifelse(df$SF_GSWLC == "Star-forming", 1, 0)
df$SF <- as.factor(df$SF)
table(df$SF)

df$LT <- ifelse(df$TType >= TType_lim, 1, 0)
df$LT <- as.factor(df$LT)
table(df$LT)

ggplot(df, aes(x = logRproj_rvir, y = logvelDisp_e, color = TType)) + 
  geom_point(size = 0.5) + 
  scale_color_distiller(palette = "RdBu", direction = 1,
                        breaks = seq(-2, 7, by = 2),
                        limits = c(min(df$TType), max(df$TType)),
                        name = "T-Type") + 
  scale_x_continuous(breaks = pretty_breaks(n = 6), name = label_logRproj_rvir) + 
  scale_y_continuous(breaks = pretty_breaks(n = 6), name = label_logvelDisp_e) + 
  theme_clean() +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 22, color = "black"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 14),
        plot.background = element_rect(color = NA)
  )

# MODELAGEM ----------------
fit_fSFG <- glm(SF ~ logvelDisp_e + logRproj_rvir, family = binomial(link = "logit"), data = df)
fit_fLTG <- glm(LT ~ logvelDisp_e + logRproj_rvir, family = binomial(link = "logit"), data = df)

summary(fit_fSFG)
summary(fit_fLTG)



summary(fit_fSFG)

coeficientes <- coef(fit_fSFG)
erro_padrao <- sqrt(diag(vcov(fit_fSFG)))  # Erro padrão calculado a partir da matriz de variância-covariância

# Calcule o z-score manualmente
z_scores <- coeficientes / erro_padrao


importance <- as.data.frame(varImp(fit_fSFG))
importance <- data.frame(overall = importance$Overall, names = rownames(importance))
importance[order(importance$overall, decreasing = T),]

importance <- as.data.frame(varImp(fit_fLTG))
importance <- data.frame(overall = importance$Overall, names = rownames(importance))
importance[order(importance$overall, decreasing = T),]

df$pred_fSFG <- predict(fit_fSFG, df, type = "response")
df$pred_fLTG <- predict(fit_fLTG, df, type = "response")

optimal_fSFG <- optimalCutoff(df$SF, df$pred_fSFG)[1]
optimal_fLTG <- optimalCutoff(df$LT, df$pred_fLTG)[1]

misClassError(df$SF, df$pred_fSFG, threshold = optimal_fSFG)
misClassError(df$LT, df$pred_fLTG, threshold = optimal_fLTG)

# -----------------------------------------------------------------------------------------

vetor_prob <- c(0.00, 1/3, 2/3, 1.00)
vetor_prob <- c(min(df$logvelDisp_e), log10(130), log10(150), max(df$logvelDisp_e))

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
      #quantil_lw_logvelDisp_e <- unname(quantile(df$logvelDisp_e, probs = prob_lw_logvelDisp_e))
      #quantil_up_logvelDisp_e <- unname(quantile(df$logvelDisp_e, probs = prob_up_logvelDisp_e))

      quantil_lw_logvelDisp_e <- prob_lw_logvelDisp_e
      quantil_up_logvelDisp_e <- prob_up_logvelDisp_e
      
      # Definindo galáxias do painel: 
      painel <- 
        df$logvelDisp_e >= quantil_lw_logvelDisp_e & 
        df$logvelDisp_e < quantil_up_logvelDisp_e 
      
      df_painel   <- df[painel,]
      ngal_painel <- nrow(df_painel)
      
      painel_label_logvelDisp_e <- paste0(
        round(10^quantil_lw_logvelDisp_e, digits = 0), " até ",
        round(10^quantil_up_logvelDisp_e, digits = 0))
      
      painel_label <- paste0("logvelDisp_e = ", painel_label_logvelDisp_e)
      
      # Modelo:
      
      # Densidades: 
      densidade_logvelDisp_e <- density(df_painel$logvelDisp_e)
      
      # Modas:
      moda_logvelDisp_e <- densidade_logvelDisp_e$x[which.max(densidade_logvelDisp_e$y)]
      
      # Medianas:
      median_logvelDisp_e <- median(df_painel$logvelDisp_e)
      
      # Médias:
      mean_logvelDisp_e <- mean(df_painel$logvelDisp_e)
      
      # ------------------------------------
      
      #bins_logRproj_rvir <- unname(quantile(df$logRproj_rvir, probs = seq(0, 1, by = 0.05)))
      bins_logRproj_rvir <- seq(min(df$logRproj_rvir), max(df$logRproj_rvir), by = 0.01)
      
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
      
      aux <- quantile(df_painel$logRproj_rvir, probs = seq(0,1, by = 0.2))
      #aux <- seq(min(df$logRproj_rvir), max(df$logRproj_rvir), by = 0.25)
      
      for (i in 1:(length(aux)-1)) {
        
        bin <- painel &
          df$logRproj_rvir >= aux[i] & 
          df$logRproj_rvir < aux[i + 1]
        
        df_bin   <- df[bin,]
        ngal_bin <- nrow(df_bin)
        
        meio_bin <- aux[i] + abs(aux[i+1] - aux[i])/2
        
        aux_label <- paste0(
          round(aux[i], digits = 3), " até ",
          round(aux[i + 1], digits = 3))
        
        ngal_SF <- sum(df_bin$SF == 1)
        ngal_Q  <- sum(df_bin$SF == 0)
        
        ngal_LT <- sum(df_bin$LT == 1)
        ngal_ET <- sum(df_bin$LT == 0)
        
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

levels <- c("50 até 130", "130 até 150","150 até 317")
levels <- c("50 até 130", "130 até 150","150 até 307")

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
