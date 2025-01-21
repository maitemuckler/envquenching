setwd("~/Work/Research/")

library(stringr)
library(InformationValue)
library(tibble)
library(data.table)
library(dplyr)
library(binom)
library(caret)
library(ggh4x)

datawd <- "Data/"
figswd <- "Figures/"

critval <- 1.96 

input_file <- "~/Work/Research/Astronomy/Data/EnvQuenching/inputModel/GSWLC/inputdata_zmax0.1_Rlim2.5_Ma12.3_flag_good==1_MANGLE_logMstar_min10.5.csv"
df <- fread(input_file)

colnames(df)[which(colnames(df) == "AGN")] <- "AGN_char"
colnames(df)[which(colnames(df) == "SF_GSWLC")] <- "SF_char"

df$AGN_char <- as.factor(df$AGN_char)
df$SF_char  <- as.factor(df$SF_char)

df$SF          <- ifelse(df$SF_char == "Star-forming", 1, 0)
df$AGN         <- ifelse(df$AGN_char == "AGN", 1, 0)

data.satellites <- subset(df, df$type == "Satellite")
data.centrals   <- subset(df, df$type == "Central")

# df <- data.satellites %>%
#   filter(AGN_char == "Non-AGN")

df <- data.satellites

rm(data.satellites, data.centrals)

# MODELAGEM

# Model 1: ----

fit_glm <- glm(SF ~ logMstar + logvelDisp_e + logRproj_rvir + TType, family = binomial(link = "logit"), data = df)

# BINNAGEM para os dados ----------------

median_bins <- data.frame(logMstar = NA,
                          logvelDisp_e = NA,
                          TType = NA,
                          
                          temp_logMstar = NA, 
                          temp_logvelDisp_e = NA,
                          temp_TType = NA,
                          temp_logRproj_rvir = NA, 
                          
                          n_total = NA,
                          n_bins = NA,
                          n_sf = NA)

model <- data.frame(logMstar = NA,
                    logvelDisp_e = NA,
                    TType = NA,
                    logRproj_rvir = NA, 
                    
                    temp_logMstar = NA, 
                    temp_logvelDisp_e = NA,
                    temp_TType = NA,
                    
                    pred = NA,
                    upr = NA,
                    lwr = NA)

logMstar_breaks      <- c(10.5, 10.7, 12.5)
logvelDisp_e_breaks  <- c(1, 2.1, 3)
TType_breaks         <- c(-3, 2.6, 7.5)
logRproj_rvir_breaks <- quantile(df$logRproj_rvir, probs = seq(0, 1, 0.2))

nbins_logMstar     <- length(logMstar_breaks)-1
nbins_logvelDisp_e <- length(logvelDisp_e_breaks)-1
nbins_TType        <- length(TType_breaks)-1

# LOOP ----------------------------
j = 1; k = 1; l = 1;

for(j in 1:nbins_logMstar){
  for(k in 1:nbins_logvelDisp_e){
    for(l in 1:nbins_TType){
    
    xx.xxm = 
      df$logMstar >= logMstar_breaks[j] & 
      df$logMstar < logMstar_breaks[j+1] &
      
      df$logvelDisp_e >= logvelDisp_e_breaks[k] & 
      df$logvelDisp_e < logvelDisp_e_breaks[k+1] &
    
      df$TType >= TType_breaks[l] & 
      df$TType < TType_breaks[l+1] 
    
    # O que vai dentro do modelo:
    temp_logMstar     <- median(df$logMstar[xx.xxm], na.rm = T)
    temp_logvelDisp_e <- median(df$logvelDisp_e[xx.xxm], na.rm = T)
    temp_TType        <- median(df$TType[xx.xxm], na.rm = T)
    logRproj          <- seq(-1, 1.35, 0.01)
    
    preds <- predict(fit_glm, data.frame(logMstar = temp_logMstar,
                                        logvelDisp_e = temp_logvelDisp_e,
                                        logRproj_rvir = logRproj,
                                        TType = temp_TType),
                    type = "link", se.fit = TRUE)
    
    logMstar_bins     <- paste0("(",format(round(logMstar_breaks[j], digits = 1), nsmall = 1), "-", format(round(logMstar_breaks[j+1], digits = 1), nsmall = 1),")")
    logvelDisp_e_bins <- paste0("(",format(round(logvelDisp_e_breaks[k], digits = 2), nsmall = 2), "-", format(round(logvelDisp_e_breaks[k+1], digits = 2), nsmall = 2),")")
    TType_bins        <- paste0("(",format(round(TType_breaks[l], digits = 2), nsmall = 2), "-", format(round(TType_breaks[l+1], digits = 2), nsmall = 2),")")
    
    model_temp <- data.frame(logMstar = logMstar_bins,
                             logvelDisp_e = logvelDisp_e_bins,
                             TType = TType_bins,
                             logRproj_rvir = logRproj,
                             
                             temp_logMstar = temp_logMstar,
                             temp_logvelDisp_e = temp_logvelDisp_e,
                             temp_TType = temp_TType,
                             
                             pred = fit_glm$family$linkinv(preds$fit),
                             upr = fit_glm$family$linkinv(preds$fit + (critval * preds$se.fit)),
                             lwr = fit_glm$family$linkinv(preds$fit - (critval * preds$se.fit)))
    
    model <- rbind(model, model_temp)
    
    for(i in 1:(length(logRproj_rvir_breaks)-1)){
      
      xx.xx = df$logRproj_rvir >= logRproj_rvir_breaks[i] & 
        df$logRproj_rvir < logRproj_rvir_breaks[i+1] & 
        xx.xxm
      
      xx.xx_SF <- df$SF[xx.xx]
      
      median_bins <- median_bins %>% 
        add_row(logMstar = logMstar_bins,
                logvelDisp_e = logvelDisp_e_bins,
                TType = TType_bins,
                
                temp_logMstar = temp_logMstar,
                temp_logvelDisp_e = temp_logvelDisp_e,
                temp_TType = temp_TType,
                temp_logRproj_rvir = median(df$logRproj_rvir[xx.xx]),
                
                n_total =  length(xx.xxm[xx.xxm== TRUE]),
                n_bins = length(xx.xx[xx.xx== TRUE]),
                n_sf = length(xx.xx_SF[xx.xx_SF == 1]))
      }
    }
  }
}

model           <- model[-1,]
median_bins     <- median_bins[-1,]
median_bins$fsf <- median_bins$n_sf/median_bins$n_bins

median_bins$confint_up <- binom.confint(x = median_bins$n_sf, n = median_bins$n_bins, 
                                        conf.level = 0.95, 
                                        method = "bayes", 
                                        type = "central")$upper
median_bins$confint_lw <- binom.confint(x = median_bins$n_sf, n = median_bins$n_bins, 
                                        conf.level = 0.95, 
                                        method = "bayes", 
                                        type = "central")$lower


# All variables as factor:
model$logMstar     <- factor(model$logMstar)
model$logvelDisp_e <- factor(model$logvelDisp_e)
model$TType        <- factor(model$TType)

median_bins$logMstar     <- factor(median_bins$logMstar)
median_bins$logvelDisp_e <- factor(median_bins$logvelDisp_e)
median_bins$TType        <- factor(median_bins$TType)

model$label       <- paste0(model$logMstar, "_", model$AGN)
median_bins$label <- paste0(median_bins$logMstar, "_", median_bins$AGN)
label_logvelDisp <- paste0("log[10]~σ~",unique(model$logvelDisp_e))

# FIGURA -------------------------------

width_figs  <- 9
height_figs <- 7
size_text_facet <- 12

ggplot() + 
  geom_line(data = model, aes(x = logRproj_rvir, y = pred, color = logMstar), linewidth = 1) + 
  geom_point(data = median_bins, aes(x = temp_logRproj_rvir, y = fsf, color = logMstar), size = 2) +
  
  geom_errorbar(data = median_bins, aes(ymin = confint_lw, ymax = confint_up, x = temp_logRproj_rvir, color = logMstar), width = 0.05, alpha = 0.4) + 
  geom_ribbon(data = model, aes(ymin = lwr, ymax = upr, x = logRproj_rvir, color = logMstar), inherit.aes = FALSE, alpha = 0.3) +
  
  facet_nested(TType  ~ logvelDisp_e) 


  #annotation_logticks(sides = "b") + 
  #scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  #coord_cartesian(xlim = 10^c(-1.28, 1.35), ylim = c(0,1)) + 

  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.position = "bottom",
        strip.text = element_text(size = size_text_facet)) + 
  
  guides(fill = guide_legend(nrow = 2, byrow = F),
         color = guide_legend(nrow = 2, byrow = F),
         linetype = guide_legend(nrow = 2, byrow = F),
         shape = guide_legend(nrow = 2, byrow = F),
         legend.spacing.y = unit(1.0, 'cm'))
