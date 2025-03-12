# Definição dos limites de dispersão de velocidade com base no redshift

zmax = 0.1

if (zmax == 0.03) {
  velDisp1 <- 50
  velDisp2 <- 75
} else if (zmax == 0.1) {
  velDisp1 <- 115
  velDisp2 <- 150
}

# Definição dos quantis
vetor_prob <- c(min(data.s$logvelDisp_e), log10(velDisp1), log10(velDisp2), max(data.s$logvelDisp_e))

# Listas para otimizar a construção dos data frames
modelos_list <- list()
bins_total_list <- list()

# Loop para calcular os modelos e bins
for (a in 1:(length(vetor_prob)-1)) {
  
  quantil_lw_logvelDisp_e <- vetor_prob[a]
  quantil_up_logvelDisp_e <- vetor_prob[a+1]
  
  painel <- data.s$logvelDisp_e >= quantil_lw_logvelDisp_e & 
    data.s$logvelDisp_e < quantil_up_logvelDisp_e 
  
  data.s_painel <- data.s[painel,]
  ngal_painel <- nrow(data.s_painel)
  
  painel_logvelDisp_e <- paste0(
    round(10^quantil_lw_logvelDisp_e, digits = 0), " até ",
    round(10^quantil_up_logvelDisp_e, digits = 0))
  
  # Criando os modelos
  bins_logRproj_rvir <- seq(min(data.s$logRproj_rvir), max(data.s$logRproj_rvir), by = 0.01)
  
  modelo <- data.frame(logvelDisp_e = median(data.s_painel$logvelDisp_e), 
                       logRproj_rvir = bins_logRproj_rvir)
  
  pred_fSFG <- predict(fit_fSFG, type = "link", se.fit = TRUE, newdata = modelo)
  pred_fLTG <- predict(fit_fLTG, type = "link", se.fit = TRUE, newdata = modelo)
  
  modelo$pred_fSFG_prob <- fit_fSFG$family$linkinv(pred_fSFG$fit)
  modelo$prob_fSFG_upr <- fit_fSFG$family$linkinv(pred_fSFG$fit + (critval * pred_fSFG$se.fit))
  modelo$prob_fSFG_lwr <- fit_fSFG$family$linkinv(pred_fSFG$fit - (critval * pred_fSFG$se.fit))
  
  modelo$pred_fLTG_prob <- fit_fLTG$family$linkinv(pred_fLTG$fit)
  modelo$prob_fLTG_upr <- fit_fLTG$family$linkinv(pred_fLTG$fit + (critval * pred_fLTG$se.fit))
  modelo$prob_fLTG_lwr <- fit_fLTG$family$linkinv(pred_fLTG$fit - (critval * pred_fLTG$se.fit))
  
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
      
      fSFG_up <- binom.confint(ngal_SF, ngal_bin, conf.level = 0.95, method = "wilson")$upper
      fSFG_lw <- binom.confint(ngal_SF, ngal_bin, conf.level = 0.95, method = "wilson")$lower
      
      fLTG_up <- binom.confint(ngal_LT, ngal_bin, conf.level = 0.95, method = "wilson")$upper
      fLTG_lw <- binom.confint(ngal_LT, ngal_bin, conf.level = 0.95, method = "wilson")$lower
      
      bins_total_list[[length(bins_total_list) + 1]] <- data.frame(
        painel_logvelDisp_e = painel_logvelDisp_e,
        meio_bin = meio_bin,
        fSFG = fSFG, fSFG_up = fSFG_up, fSFG_lw = fSFG_lw,
        fLTG = fLTG, fLTG_up = fLTG_up, fLTG_lw = fLTG_lw)
    }
  }
}

# Convertendo listas para data frames
modelos <- do.call(rbind, modelos_list)
bins_total <- do.call(rbind, bins_total_list)

# Garantindo que painel_logvelDisp_e é um fator
modelos$painel_logvelDisp_e <- factor(modelos$painel_logvelDisp_e)
bins_total$painel_logvelDisp_e <- factor(bins_total$painel_logvelDisp_e)

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
  
  facet_grid(. ~ painel_logvelDisp_e, scales = "free_x") + 
  
  ylab("Fraction") + 
  xlab(label_logRproj_rvir) + 
  
  scale_color_manual(values = c("fSFG" = "#4361ee", "fLTG" = "#e63946")) + 
  scale_fill_manual(values = c("fSFG" = "#4361ee", "fLTG" = "#e63946")) + 
  
  theme_bw() +
  theme(text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = "bottom")
