setwd("~/Work/Research/")

library(stringr)
library(InformationValue)
library(tibble)
library(data.table)
library(dplyr)
library(binom)
library(caret)

datawd <- "Data/"
figswd <- "Figures/"

critval <- 1.96 

# Ma_vec   <- c(12.3, 13)
# cat_vec  <- c("GSWLC", "MPA-JHU")
# zmax_vec <- c(0.03, 0.1)
# lim_vec  <- c("_nocentral", "")

Ma_vec   <- c(12.3)
cat_vec  <- c("MPA-JHU")
zmax_vec <- c(0.1)
lim_vec  <- c("")

for (lim_op in lim_vec) {
  for (catalogue_op in cat_vec) {
    for (z_max in zmax_vec) {
      for (Mhalo_min_group in Ma_vec) {

        if(catalogue_op == "GSWLC"){SF_name = "SF_gswlc"}
        if(catalogue_op == "MPA-JHU"){SF_name = "SF_mpajhu"}
        
        input_file <- paste0(datawd,
                             "input_model/",
                             "input_model_",
                             catalogue_op,
                             "_zmax", z_max,
                             "_Ma", Mhalo_min_group,
                             "_flag_good==1_MANGLE",
                             lim_op,
                             ".csv")
        
        gal <- read.csv(input_file)
        
        if(z_max == 0.03){gal <- subset(gal, gal$logMstar < 10.5)}
        
        colnames(gal)[which(colnames(gal) == "AGN")] <- "AGN_char"
        colnames(gal)[which(colnames(gal) == SF_name)] <- "SF_char"
        
        gal$SF          <- ifelse(gal$SF_char == "Star forming", 1, 0)
        gal$AGN         <- ifelse(gal$AGN_char == "AGN", 1, 0)
        
        data.model    <- subset(gal, gal$type == "Satellite")
        data.centrais <- subset(gal, gal$type == "Central")
        
        fit_glm <- glm(SF ~ 
                         poly(logMgroup,     1, raw = T) +
                         poly(logMstar,      3, raw = T) +
                         poly(logvelDisp_e,  3, raw = T) +
                         poly(logRproj_rvir, 1, raw = T) + 
                         AGN,
                       
                       family = binomial(link = "logit"), 
                       data = data.model)
        
        # BINNAGEM para os dados ----------------
        
        median_bins <- data.frame(logMgroup = NA,
                                  logMstar = NA,
                                  logvelDisp_e = NA,
                                  AGN = NA,
                                  
                                  temp_mhalo = NA, 
                                  temp_mstar = NA, 
                                  temp_sigma_re = NA,
                                  temp_rbins = NA, 
                                  
                                  n_total = NA,
                                  n_bins = NA,
                                  n_sf = NA,
                                  
                                  n_centrais = NA,
                                  n_sf_c = NA)
        
        model <- data.frame(logMgroup = NA,
                            logMstar = NA,
                            logRproj = NA, 
                            logvelDisp_e = NA,
                            AGN = NA,
                            temp_mhalo = NA,
                            temp_mstar = NA,
                            temp_sigma_re = NA,
                            pred = NA,
                            upr = NA,
                            lwr = NA)
        
        if(z_max == 0.03){mass_breaks = quantile(data.model$logMstar, probs = c(0, 1))}
        if(z_max == 0.1){mass_breaks = c(10.5, 10.67, 12.5)}
        
        if(z_max == 0.03){mhalo_breaks = quantile(data.model$logMgroup, probs = c(0, 1))}
        if(z_max == 0.1 & Mhalo_min_group == 12.3){mhalo_breaks = c(12.3, 14, 15.5)}
        if(z_max == 0.1 & Mhalo_min_group == 13){mhalo_breaks = c(13, 14, 15.5)}
        
        if(z_max == 0.03){velDisp_breaks = c(1, 2.1, 3)}
        if(z_max == 0.1){velDisp_breaks = c(1, 1.9, 2.1, 3)}
        
        rbins          <- c(0, quantile(data.model$Rproj_rvir, seq(0, 1, by = 0.2)))
        
        nbins_mhalo    <- length(mhalo_breaks)-1
        nbins_mgal     <- length(mass_breaks)-1
        nbins_sigma_re <- length(velDisp_breaks)-1
        agn_bins       <- c(0,1)
        
        # i=1
        # j=1
        # l=1
        # g=0
        
        for (g in agn_bins) {
          for(i in 1:nbins_mhalo){
            for(j in 1:nbins_mgal){
              for(l in 1:nbins_sigma_re){
                
                xx.xxm = 
                  data.model$logMstar >= mass_breaks[j] & 
                  data.model$logMstar < mass_breaks[j+1] &
                  
                  data.model$logMgroup >= mhalo_breaks[i] & 
                  data.model$logMgroup < mhalo_breaks[i+1] &
                  
                  data.model$logvelDisp_e >= velDisp_breaks[l] & 
                  data.model$logvelDisp_e < velDisp_breaks[l+1] &
                  
                  data.model$AGN == g
                
                # CENTRAIS
                xx.xxm_c = 
                  data.centrais$logMstar >= mass_breaks[j] & 
                  data.centrais$logMstar < mass_breaks[j+1] &
                  
                  data.centrais$logMgroup >= mhalo_breaks[i] & 
                  data.centrais$logMgroup < mhalo_breaks[i+1] &
                  
                  data.centrais$logvelDisp_e >= velDisp_breaks[l] & 
                  data.centrais$logvelDisp_e < velDisp_breaks[l+1] &
                  data.centrais$AGN == g
                
                
                temp_mhalo    = median(data.model$logMgroup[xx.xxm], na.rm = T)
                temp_mstar    = median(data.model$logMstar[xx.xxm], na.rm = T)
                temp_sigma_re = median(data.model$logvelDisp_e[xx.xxm], na.rm = T)
                
                logRproj = seq(-1, 1.35, 0.01)
                
                preds = predict(fit_glm, data.frame(logMgroup = temp_mhalo,
                                                    logMstar = temp_mstar,
                                                    logRproj_rvir = logRproj,
                                                    logvelDisp_e = temp_sigma_re,
                                                    AGN = g),
                                type = "link", se.fit = TRUE)
                
                logMgroup_bins    <- paste0("(",format(round(mhalo_breaks[i], digits = 1), nsmall = 1), "-",
                                            format(round(mhalo_breaks[i+1], digits = 1), nsmall = 1),")")
                logMstar_bins     <- paste0("(",format(round(mass_breaks[j], digits = 1), nsmall = 1), "-",
                                            format(round(mass_breaks[j+1], digits = 1), nsmall = 1),")")
                logvelDisp_e_bins <- paste0("(",format(round(velDisp_breaks[l], digits = 2), nsmall = 2), "-",
                                            format(round(velDisp_breaks[l+1], digits = 2), nsmall = 2),")")
                
                model_temp <- data.frame(logMgroup = logMgroup_bins,
                                         logMstar = logMstar_bins,
                                         logRproj,
                                         logvelDisp_e = logvelDisp_e_bins,
                                         AGN = g,
                                         temp_mhalo = temp_mhalo,
                                         temp_mstar = temp_mstar,
                                         temp_sigma_re = temp_sigma_re,
                                         pred = fit_glm$family$linkinv(preds$fit),
                                         upr = fit_glm$family$linkinv(preds$fit + (critval * preds$se.fit)),
                                         lwr = fit_glm$family$linkinv(preds$fit - (critval * preds$se.fit)))
                
                model <- rbind(model, model_temp)
                
                for(k in 1:(length(rbins)-1)){
                  
                  xx.xx = data.model$Rproj_rvir >= rbins[k] & 
                    data.model$Rproj_rvir < rbins[k+1] & 
                    xx.xxm
                  
                  xx.xx_SF <- data.model$SF[xx.xx]
                  
                  # CENTRAIS
                  xx.xx_SF_c <- data.centrais$SF[xx.xxm_c]
                  
                  median_bins <- median_bins %>% 
                    add_row(logMgroup = logMgroup_bins,
                            logMstar = logMstar_bins,
                            logvelDisp_e = logvelDisp_e_bins,
                            AGN = g,
                            temp_mhalo = temp_mhalo,
                            temp_mstar = temp_mstar,
                            temp_sigma_re = temp_sigma_re,
                            temp_rbins = median(data.model$Rproj_rvir[xx.xx]),
                            n_total =  length(xx.xxm[xx.xxm== TRUE]),
                            n_bins = length(xx.xx[xx.xx== TRUE]),
                            n_sf = length(xx.xx_SF[xx.xx_SF == 1]),
                            n_centrais = length(xx.xxm_c[xx.xxm_c== TRUE]),
                            n_sf_c = length(xx.xx_SF_c[xx.xx_SF_c == 1]))
                  
                }
              }
            }
          }
        }
        
        median_bins <- median_bins[-1,]
        nrow(median_bins)/3
        median_bins$fsf <- median_bins$n_sf/median_bins$n_bins
        median_bins <- median_bins %>%
          na.omit()
        
        model <- na.omit(model) # remover a primeira linha que é NA
        
        # centrais
        median_bins$fsf_c <- median_bins$n_sf_c/median_bins$n_centrais
        
        if(any(is.na(median_bins$fsf_c))){
          median_bins$fsf_c[which(is.na(median_bins$fsf_c))] <- 0
        }
        
        median_bins$confint_up <- binom.confint(x = median_bins$n_sf, n = median_bins$n_bins, 
                                                conf.level = 0.95, 
                                                method = "bayes", 
                                                type = "central")$upper
        median_bins$confint_lw <- binom.confint(x = median_bins$n_sf, n = median_bins$n_bins, 
                                                conf.level = 0.95, 
                                                method = "bayes", 
                                                type = "central")$lower
        
        # CENTRAIS
        median_bins$confint_up_c <- binom.confint(x = median_bins$n_sf_c, n = median_bins$n_centrais, 
                                                  conf.level = 0.95, 
                                                  method = "bayes", 
                                                  type = "central")$upper
        median_bins$confint_lw_c <- binom.confint(x = median_bins$n_sf_c, n = median_bins$n_centrais, 
                                                  conf.level = 0.95, 
                                                  method = "bayes", 
                                                  type = "central")$lower
        
        label_logvelDisp <- paste0("log[10]~σ~",unique(model$logvelDisp_e))
        label_logMgroup  <- paste0("log[10]~M['h']~",unique(model$logMgroup))
        label_whan       <- c("Non-AGN", "AGN")
        
        # All variables as factor:
        model$logMstar       <- factor(model$logMstar)
        model$logvelDisp_e   <- factor(model$logvelDisp_e, labels = label_logvelDisp)
        model$logMgroup      <- factor(model$logMgroup, labels = label_logMgroup)
        model$AGN       <- factor(model$AGN, labels = label_whan)
        
        median_bins$logMstar     <- factor(median_bins$logMstar)
        median_bins$logvelDisp_e <- factor(median_bins$logvelDisp_e, labels = label_logvelDisp)
        median_bins$logMgroup    <- factor(median_bins$logMgroup, labels = label_logMgroup)
        median_bins$AGN     <- factor(median_bins$AGN, labels = label_whan)
        
        model$label       <- paste0(model$logMstar, "_", model$AGN)
        median_bins$label <- paste0(median_bins$logMstar, "_", median_bins$AGN)
    
        
        write.csv(median_bins, row.names = F,
                  file = paste0(datawd,
                                "output_model/median_bins/",
                                "median_bins_",
                                catalogue_op,
                                "_zmax", z_max,
                                "_Ma", Mhalo_min_group,
                                "_flag_good==1_MANGLE",
                                lim_op,
                                ".csv"))
        
        write.csv(model, row.names = F,
                  file = paste0(datawd,
                                "output_model/model/",
                                "model_",
                                catalogue_op,
                                "_zmax", z_max,
                                "_Ma", Mhalo_min_group,
                                "_flag_good==1_MANGLE",                                
                                lim_op,
                                ".csv"))

      }
    }
  }
}



