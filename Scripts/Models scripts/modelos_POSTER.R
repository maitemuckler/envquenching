## Carregando pacotes ----
library(stringr)
library(data.table)
library(dplyr)
library(binom)
library(InformationValue)
library(caret)

## Diretórios ----
wdmain       <- "~/Work/Research/Astronomy/Projects/environmental-quenching/"
wdcode       <- "~/Work/Research/Astronomy/Projects/environmental-quenching/Scripts/"
wdinputdata  <- "~/Work/Research/Astronomy/Data/EnvQuenching/inputModel/"
wdoutputdata <- "~/Work/Research/Astronomy/Data/EnvQuenching/outputModel/"
wdimportance <- "~/Work/Research/Astronomy/Data/EnvQuenching/importanceModel/"
wdfigs       <- "~/Work/Research/Astronomy/Projects/environmental-quenching/Figures/"

## Seed ----
set.seed(123)

## Definindo opções ----
catalog <- "GSWLC"
critval <- 1.96 
zmax    <- 0.1
Rlim    <- 2.5
Ma      <- 12.3

tipo_de_modelo <- "NAGN"

## Definindo input e output files ----
input_file             <- paste0("inputdata_zmax", zmax, "_Rlim", Rlim, "_Ma", Ma, "_flag_good==1_MANGLE_logMstar_min10.5.csv")
output_file_bins       <- paste0("bins_", catalog, "_zmax", zmax, "_Ma", Ma, ".csv")
output_file_model      <- paste0("model_", catalog, "_zmax", zmax, "_Ma", Ma, ".csv")
output_file_importance <- paste0("importance_", catalog, "_zmax", zmax, "_Ma", Ma, ".csv")

## Lendo os dados ----
df  <- fread(paste0(wdinputdata, catalog, "/", input_file)) 

## Definindo variáveis de acordo com catálogo ----
if(catalog == "GSWLC"){
  Mstar_name <- "logMstar"
  SF_name    <- "SF_GSWLC"
}

if(catalog == "MPA-JHU"){
  Mstar_name <- "lgm_tot_p50"
  SF_name    <- "SF_MPAJHU"
}

## Selecionando variáveis ----
variaveis <- c("type",
               SF_name,
               Mstar_name,
               "AGN",
               "logvelDisp_e",
               "logMgroup",
               "Rproj_rvir",
               "logRproj_rvir")

df <- df %>%
  select(all_of(variaveis))

## Renomear variáveis ----
df <- df %>%
  rename(SF = SF_name) %>%
  rename(logMstar = Mstar_name)

## Definindo variáveis factor como 0 e 1 ----
df$SF  <- ifelse(df$SF == "Star-forming", 1, 0)
df$SF  <- as.factor(df$SF)

## Separando dados em satélites e centrais ----
data.satellites <- subset(df, df$type == "Satellite")
data.centrals   <- subset(df, df$type == "Central")

## Separando dados em AGN e não-AGN -----
s_agn  <- subset(data.satellites, data.satellites$AGN == "AGN")
s_nagn <- subset(data.satellites, data.satellites$AGN == "Non-AGN")

c_agn  <- subset(data.centrals, data.centrals$AGN == "AGN")
c_nagn <- subset(data.centrals, data.centrals$AGN == "Non-AGN")

## Removendo variáveis e data.frames desnecessários ----
rm(df, data.satellites, data.centrals)

s_agn  <- s_agn %>% select(-c(type, AGN))
s_nagn <- s_nagn %>% select(-c(type, AGN))
c_agn  <- c_agn %>% select(-c(type, AGN))
c_nagn <- c_nagn %>% select(-c(type, AGN))

## Modelos ----

if(tipo_de_modelo == "AGN"){
  df     <- s_agn
  df_c   <- c_agn
}

if(tipo_de_modelo == "NAGN"){
  df     <- s_nagn
  df_c   <- c_nagn
}

modelo <- glm(SF ~ logMstar + logvelDisp_e + logRproj_rvir + logMgroup, family = binomial(link = "logit"), data = df)

summary(modelo)

importance <- as.data.frame(varImp(modelo))
importance <- data.frame(overall = importance$Overall, names = rownames(importance))
importance[order(importance$overall, decreasing = T),]

write.csv(importance, paste0(wdimportance, catalog, "/", tipo_de_modelo, "/", output_file_importance), row.names = F)

df$pred <- predict(modelo, df, type = "response")
#optimal    <- optimalCutoff(s_agn$SF, s_agn$pred)[1]
optimal <- 0.5

misClassError(df$SF, df$pred, threshold = optimal)

## Binnagem ----

median_bins <- data.frame(logMgroup = NA,
                          logMstar = NA,
                          logvelDisp_e = NA,
                          
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
                    temp_mhalo = NA,
                    temp_mstar = NA,
                    temp_sigma_re = NA,
                    pred = NA,
                    upr = NA,
                    lwr = NA)

if(zmax == 0.03){mass_breaks = quantile(df$logMstar, probs = c(0, 1))}
if(zmax == 0.1){mass_breaks = c(10.5, 10.67, 12.5)}

if(zmax == 0.03){mhalo_breaks = quantile(df$logMgroup, probs = c(0, 1))}
if(zmax == 0.1 & Ma == 12.3){mhalo_breaks = c(12.3, 14, 15.5)}
if(zmax == 0.1 & Ma == 13){mhalo_breaks = c(13, 14, 15.5)}

if(zmax == 0.03){velDisp_breaks = c(1, 2.1, 3)}
if(zmax == 0.1){velDisp_breaks = c(1, 1.9, 2.1, 3)}

rbins          <- c(0, quantile(df$Rproj_rvir, seq(0, 1, by = 0.2)))

nbins_mhalo    <- length(mhalo_breaks)-1
nbins_mgal     <- length(mass_breaks)-1
nbins_sigma_re <- length(velDisp_breaks)-1

i=1
j=1
l=1

for(i in 1:nbins_mhalo){
  for(j in 1:nbins_mgal){
    for(l in 1:nbins_sigma_re){
      
      xx.xxm = 
        df$logMstar >= mass_breaks[j] & 
        df$logMstar < mass_breaks[j+1] &
        
        df$logMgroup >= mhalo_breaks[i] & 
        df$logMgroup < mhalo_breaks[i+1] &
        
        df$logvelDisp_e >= velDisp_breaks[l] & 
        df$logvelDisp_e < velDisp_breaks[l+1] 
      
      # CENTRAIS
      xx.xxm_c = 
        df_c$logMstar >= mass_breaks[j] & 
        df_c$logMstar < mass_breaks[j+1] &
        
        df_c$logMgroup >= mhalo_breaks[i] & 
        df_c$logMgroup < mhalo_breaks[i+1] &
        
        df_c$logvelDisp_e >= velDisp_breaks[l] & 
        df_c$logvelDisp_e < velDisp_breaks[l+1] 
      
      
      temp_mhalo    = median(df$logMgroup[xx.xxm], na.rm = T)
      temp_mstar    = median(df$logMstar[xx.xxm], na.rm = T)
      temp_sigma_re = median(df$logvelDisp_e[xx.xxm], na.rm = T)
      
      logRproj = seq(-1, 1.35, 0.01)
      
      preds = predict(modelo, data.frame(logMgroup = temp_mhalo,
                                         logMstar = temp_mstar,
                                         logRproj_rvir = logRproj,
                                         logvelDisp_e = temp_sigma_re),
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
                               temp_mhalo = temp_mhalo,
                               temp_mstar = temp_mstar,
                               temp_sigma_re = temp_sigma_re,
                               pred = modelo$family$linkinv(preds$fit),
                               upr = modelo$family$linkinv(preds$fit + (critval * preds$se.fit)),
                               lwr = modelo$family$linkinv(preds$fit - (critval * preds$se.fit)))
      
      model <- rbind(model, model_temp)
      
      for(k in 1:(length(rbins)-1)){
        
        xx.xx = df$Rproj_rvir >= rbins[k] & 
          df$Rproj_rvir < rbins[k+1] & 
          xx.xxm
        
        xx.xx_SF <- df$SF[xx.xx]
        
        # CENTRAIS
        xx.xx_SF_c <- df_c$SF[xx.xxm_c]
        
        median_bins <- median_bins %>% 
          add_row(logMgroup = logMgroup_bins,
                  logMstar = logMstar_bins,
                  logvelDisp_e = logvelDisp_e_bins,
                  temp_mhalo = temp_mhalo,
                  temp_mstar = temp_mstar,
                  temp_sigma_re = temp_sigma_re,
                  temp_rbins = median(df$Rproj_rvir[xx.xx]),
                  n_total =  length(xx.xxm[xx.xxm== TRUE]),
                  n_bins = length(xx.xx[xx.xx== TRUE]),
                  n_sf = length(xx.xx_SF[xx.xx_SF == 1]),
                  n_centrais = length(xx.xxm_c[xx.xxm_c== TRUE]),
                  n_sf_c = length(xx.xx_SF_c[xx.xx_SF_c == 1]))
        
      }
    }
  }
}

## Arumando data.frame dos bins ----
median_bins <- median_bins %>%
  filter(!if_all(everything(), is.na))

nrow(median_bins)/3

median_bins$fsf   <- median_bins$n_sf/median_bins$n_bins
median_bins$fsf_c <- median_bins$n_sf_c/median_bins$n_centrais # centrais
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

## Arrumando data.frame do modelo ----
model <- model %>%
  filter(!if_all(everything(), is.na))

## Arrumando labels ----
label_logvelDisp <- paste0("log[10]~σ~",unique(model$logvelDisp_e))
label_logMgroup  <- paste0("log[10]~M['h']~",unique(model$logMgroup))

## Definindo todas variáveis como factor ----
model$logMstar       <- factor(model$logMstar)
model$logvelDisp_e   <- factor(model$logvelDisp_e, labels = label_logvelDisp)
model$logMgroup      <- factor(model$logMgroup, labels = label_logMgroup)

median_bins$logMstar     <- factor(median_bins$logMstar)
median_bins$logvelDisp_e <- factor(median_bins$logvelDisp_e, labels = label_logvelDisp)
median_bins$logMgroup    <- factor(median_bins$logMgroup, labels = label_logMgroup)

## Salvando tabela ----
write.csv(median_bins, file = paste0(wdoutputdata, catalog, "/", tipo_de_modelo, "/", output_file_bins), row.names = F)
write.csv(model, file = paste0(wdoutputdata, catalog, "/", tipo_de_modelo, "/", output_file_model), row.names = F)
