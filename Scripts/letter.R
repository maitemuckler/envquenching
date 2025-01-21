setwd("~/Work/Research/Astronomy/Projects/environmental-quenching/")

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

wddata <- "~/Work/Research/Astronomy/Data/environmental-quenching-data/"

zmax          <- 0.1
catalog       <- "GSWLC"
TType_lim     <- 0
critval       <- 1.96 
mass_complete <- "y"
Ma            <- 12.3

SF_name    <- ifelse(catalog == "GSWLC", "SF_GSWLC", "SF_MPAJHU")

if(mass_complete == "y"){
  input_data <- paste0("/inputdata_zmax", zmax, "_Rlim2.5_Ma", Ma, "_flag_good==1_MANGLE_logMstar_min10.5.csv")
}

if(mass_complete == "n"){
  input_data <- paste0("/inputdata_zmax", zmax, "_Rlim2.5_Ma", Ma, "_flag_good==1_MANGLE.csv")
}

df    <- fread(paste0(wddata, "inputModel/", catalog, input_data))

if(mass_complete == "n"){ # fazer Vmax == 0 ser 1e-6 e normalizar os pesos
  df$Vmax[df$Vmax == 0] <- 1e-6
  df$weights <- 1/df$Vmax
  df$weights <- df$weights / max(df$weights)
}

colnames(df)[which(colnames(df) == "AGN")]   <- "AGN_char"
colnames(df)[which(colnames(df) == SF_name)] <- "SF_char"

df$AGN_char <- as.factor(df$AGN_char)
df$SF_char  <- as.factor(df$SF_char)

df$SF          <- ifelse(df$SF_char == "Star-forming", 1, 0)
df$AGN         <- ifelse(df$AGN_char == "AGN", 1, 0)

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
ggplot(data_binned.s, aes(x = median_logRproj_rvir, y = median_TType, color = factor(SF))) +
  geom_smooth(aes(group = SF))          +
  facet_grid(. ~ sigma_bin) + 
  labs(x = "logRproj", y = "t-type", color = "Class") +
  scale_color_manual(values = c("red", "blue")) 

ggplot(data_binned.s, aes(x = median_logRproj_rvir, y = median_TType, color = factor(SF))) +
  geom_smooth(aes(group = SF)) +
  labs(x = "logRproj", y = "t-type", color = "Class") +
  scale_color_manual(values = c("red", "blue")) 



# Modelos ------------------------------------------------------------------------------------------------

if(mass_complete == "y"){
  fit_fSFG <- glm(SF ~ logvelDisp_e + logRproj_rvir, family = binomial(link = "logit"), data = data.s)
  fit_fLTG <- glm(LT ~ logvelDisp_e + logRproj_rvir, family = binomial(link = "logit"), data = data.s)
}

if(mass_complete == "n"){
  fit_fSFG <- glm(SF ~ logvelDisp_e + logRproj_rvir, family = binomial(link = "logit"), data = data.s, weights = weights)
  fit_fLTG <- glm(LT ~ logvelDisp_e + logRproj_rvir, family = binomial(link = "logit"), data = data.s, weights = weights)
}

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
