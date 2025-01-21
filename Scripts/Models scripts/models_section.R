setwd("~/Work/Research/")

library(stringr)
library(InformationValue)
library(tibble)
library(data.table)
library(dplyr)
library(binom)
library(caret)
library(pROC)

datawd <- "Data/"
figswd <- "Figures/"

critval <- 1.96 

input_file <- "~/Work/Research/Astronomy/Data/InputModel/inputdata_zmax0.03_Rlim2.5_Ma12.3_flag_good==1_MANGLE_logMstar_min10.5.csv"
df <- fread(input_file)

colnames(df)[which(colnames(df) == "AGN")] <- "AGN_char"
colnames(df)[which(colnames(df) == "SF_GSWLC")] <- "SF_char"

df$AGN_char <- as.factor(df$AGN_char)
df$SF_char  <- as.factor(df$SF_char)

df$SF          <- ifelse(df$SF_char == "Star-forming", 1, 0)
df$AGN         <- ifelse(df$AGN_char == "AGN", 1, 0)

data.satellites <- subset(df, df$type == "Satellite")
data.centrals   <- subset(df, df$type == "Central")

#save(data.satellites, file = "~/Work/Research/Astronomy/Projects/environmental-quenching/Data/satelites.RData", compress = "xz")
#save(data.centrals, file = "~/Work/Research/Astronomy/Projects/environmental-quenching/Data/centrais.RData", compress = "xz")

# Binnagem:
data.satellites <- data.satellites %>% 
  mutate(logMstar_BIN            = cut_number(logMstar,           n = 5, right = F)) %>%
  mutate(Rproj_rvir_BIN          = cut_number(Rproj_rvir,         n = 5, right = F)) %>%
  mutate(logMgroup_BIN           = cut_number(logMgroup,          n = 5, right = F)) %>%
  mutate(logvelDisp_e_BIN        = cut_number(logvelDisp_e,       n = 5, right = F)) %>%
  mutate(vlos_vvir_BIN           = cut_number(vlos_vvir,          n = 5, right = F)) %>%
  mutate(d4000_n_BIN             = cut_number(d4000_n,            n = 5, right = F)) %>%
  # AGN
  
  mutate(logdens_n_Rproj0.5_BIN  = cut_number(logdens_n_Rproj0.5, n = 3, right = F)) %>%
  mutate(logdens_n_Rproj1_BIN    = cut_number(logdens_n_Rproj1,   n = 3, right = F)) %>%
  mutate(Ngals_vol_Rproj0.5_BIN  = cut_number(Ngals_vol_Rproj0.5, n = 3, right = F)) %>%
  mutate(Ngals_vol_Rproj1_BIN    = cut_number(Ngals_vol_Rproj1,   n = 3, right = F)) %>%
  mutate(logdens_proj_Neq1_BIN   = cut_number(logdens_proj_Neq1,  n = 3, right = F)) %>%
  mutate(logdens_proj_Neq3_BIN   = cut_number(logdens_proj_Neq3,  n = 3, right = F)) %>%
  mutate(logdens_proj_Neq5_BIN   = cut_number(logdens_proj_Neq5,  n = 3, right = F)) %>%
  
  mutate(nb_BIN            = cut_number(nb,           n = 5, right = F)) %>%
  mutate(B_T_r_BIN         = cut_number(B_T_r,        n = 5, right = F)) %>%
  mutate(e_BIN             = cut_number(e,            n = 5, right = F)) %>%
  mutate(P_disk_BIN        = cut_number(P_disk,       n = 5, right = F)) %>%
  mutate(P_edge_on_BIN     = cut_number(P_edge_on,    n = 5, right = F)) %>%
  mutate(P_bar_GZ2_BIN     = cut_number(P_bar_GZ2,    n = 5, right = F)) %>%
  mutate(P_bar_Nair10_BIN  = cut_number(P_bar_Nair10, n = 5, right = F)) %>%
  mutate(P_merg_BIN        = cut_number(P_merg,       n = 5, right = F)) %>%
  mutate(P_bulge_BIN       = cut_number(P_bulge,      n = 5, right = F)) %>%
  mutate(P_cigar_BIN       = cut_number(P_cigar,      n = 5, right = F)) %>%
  mutate(P_S0_BIN          = cut_number(P_S0,         n = 5, right = F)) %>%
  
  mutate(distLine_GSWLC_BIN = cut_number(distLine_GSWLC, n = 5, right = F)) %>%
  mutate(logSigma_SFR_BIN   = cut_number(logSigma_SFR,   n = 5, right = F)) %>%
  
  mutate(logMstar_central_BIN        = cut_number(logMstar_central,       n = 5, right = F)) %>%
  mutate(logvelDisp_e_central_BIN    = cut_number(logvelDisp_e_central,   n = 5, right = F)) %>%
  mutate(vlos_vvir_central_BIN       = cut_number(vlos_vvir_central,      n = 5, right = F)) %>%
  mutate(distLine_GSWLC_central_BIN  = cut_number(distLine_GSWLC_central, n = 5, right = F)) %>%
  mutate(P_disk_central_BIN          = cut_number(P_disk_central,         n = 5, right = F)) %>%
  mutate(P_edge_on_central_BIN       = cut_number(P_edge_on_central,      n = 5, right = F)) %>%
  mutate(P_bar_GZ2_central_BIN       = cut_number(P_bar_GZ2_central,      n = 5, right = F)) %>%
  mutate(P_bar_Nair10_central_BIN    = cut_number(P_bar_Nair10_central,   n = 5, right = F)) %>%
  mutate(P_merg_central_BIN          = cut_number(P_merg_central,         n = 5, right = F)) %>%
  mutate(P_bulge_central_BIN         = cut_number(P_bulge_central,        n = 5, right = F)) %>%
  mutate(B_T_r_central_BIN           = cut_number(B_T_r_central,          n = 5, right = F)) %>%
  mutate(e_central_BIN               = cut_number(e_central,              n = 5, right = F)) %>%
  mutate(P_cigar_central_BIN         = cut_number(P_cigar_central,        n = 5, right = F)) %>%
  mutate(P_bulge_central_BIN         = cut_number(P_bulge_central,        n = 5, right = F))

# Calcular a mediana de cada intervalo
medianas <- data.satellites %>%
  group_by(Rproj_rvir_BIN) %>%
  summarise(mediana = median(Rproj_rvir))

# Criar os rótulos baseados nas medianas
labels_medianas <- as.character(medianas$mediana)

# Usar cut_number novamente, mas agora definindo os rótulos como as medianas
data.satellites <- data.satellites %>%
  mutate(Rproj_rvir_BIN_med = cut_number(Rproj_rvir, n = 5, right = F, labels = labels_medianas))

trainIndex <- createDataPartition(data.satellites$SF, p = 0.7, list = FALSE)
train_data <- data.satellites[trainIndex, ]
test_data  <- data.satellites[-trainIndex, ]

rm(data.satellites)
rm(data.centrals)
rm(trainIndex)
rm(df)

# Model 1: ----

model1 <- glm(SF ~ TType + logvelDisp_e + conc + logRproj_rvir + logdens_proj_Neq1, family = binomial(link = "logit"), data = train_data)
summary(model1)

test_data$prob1    <- predict(model1, newdata = test_data, type = "response")
optimal1           <- optimalCutoff(test_data$SF, test_data$prob1)[1]
test_data$SF_pred1 <- ifelse(test_data$prob1 > optimal1, 1, 0) 

test_data$SF_pred1 <- factor(test_data$SF_pred1)
test_data$SF       <- factor(test_data$SF)

confusion_matrix <- confusionMatrix(test_data$SF_pred1, test_data$SF)
confusion_matrix$overall['Accuracy']

misClassError(test_data$SF, test_data$prob1, threshold = optimal1) 

roc_obj1 <- roc(test_data$SF, test_data$prob1) 
plot(roc_obj1, main = "Curva ROC", col = "blue")

imp1 <- as.data.frame(varImp(model1))
imp1 <- data.frame(overall = imp1$Overall, names = rownames(imp1))
imp1[order(imp1$overall, decreasing = T),]

bins <- test_data %>%
  group_by(logMstar_BIN, logvelDisp_e_BIN, Rproj_rvir_BIN_med, SF_char) %>% 
  summarise(n = n()) %>%
  mutate(fraçao = n / sum(n) * 100) %>%
  ungroup()

ggplot(test_data, aes(x = Rproj_rvir_BIN_med, y = prob1)) + 
  geom_line() + 
  facet_wrap(logMstar_BIN ~ logvelDisp_e_BIN)






























model1 <- glm(SF ~ logvelDisp_e, family = binomial(link = "logit"), data = train_data)
summary(model1)




# Model 2: ----



test_data$prob2    <- predict(model2, newdata = test_data, type = "response")
optimal2           <- optimalCutoff(test_data$SF, test_data$prob2)[1] 
test_data$SF_pred2 <- ifelse(test_data$prob2 > optimal2, 1, 0) 
confusionMatrix(data = as.factor(test_data$SF_pred2), reference = as.factor(test_data$SF)) 
misClassError(test_data$SF, test_data$prob2, threshold = optimal2) 

roc_obj2 <- roc(test_data$SF, test_data$prob2) 
plot(roc_obj2, main = "Curva ROC", col = "blue")

imp2 <- as.data.frame(varImp(model2))
imp2 <- data.frame(overall = imp2$Overall, names = rownames(imp2))
imp2[order(imp2$overall, decreasing = T),]











# Configurar controle de treino para validação cruzada de 10 folds
train_control <- trainControl(method = "cv", number = 10)

# Ajustar o modelo usando validação cruzada
modelo_cv <- train(SF ~ logvelDisp_e_scale, data = train_data, method = "glm", family = "binomial", trControl = train_control)

# Resultados da validação cruzada
print(modelo_cv)

summary(model1)$coefficients












# Gráficos ----




# Model 2:

test_data <- test_data %>% 
  mutate(logMstar_BIN   = cut_number(logMstar,     n = 2, right = F)) %>%
  mutate(Rproj_rvir_BIN = cut_number(Rproj_rvir,   n = 5, right = F))

# sigma
data.satellites$logvelDisp_BIN <- NA
data.satellites$logvelDisp_BIN[which(data.satellites$logvelDisp_e < 1.9)] <- "[1.0,1.9)"
data.satellites$logvelDisp_BIN[which(data.satellites$logvelDisp_e >= 1.9 & 
                                       data.satellites$logvelDisp_e < 2.1)] <- "[1.9,2.1)"
data.satellites$logvelDisp_BIN[which(data.satellites$logvelDisp_e >= 2.1)] <- "[2.1,3.0)"
data.satellites$logvelDisp_BIN <- as.factor(data.satellites$logvelDisp_BIN)

# mhalo
data.satellites$logMgroup_BIN <- NA
data.satellites$logMgroup_BIN[which(data.satellites$logMgroup < 13)] <- "[12.3,13)"
data.satellites$logMgroup_BIN[which(data.satellites$logMgroup >= 13 & 
                                      data.satellites$logMgroup < 14)] <- "[13,14)"
data.satellites$logMgroup_BIN[which(data.satellites$logMgroup >= 14)] <- "[14,15.5)"
data.satellites$logMgroup_BIN <- as.factor(data.satellites$logMgroup_BIN)