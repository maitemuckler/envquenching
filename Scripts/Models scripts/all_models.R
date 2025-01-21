library(data.table)
library(caret)
library(InformationValue)
library(dplyr)

### Definir diretórios ----
wdmain       <- "~/Work/Research/Astronomy/"
wdproject    <- paste0(wdmain, "Projects/environmental-quenching/")
wdcode       <- paste0(wdproject, "Scripts/")
wddata       <- paste0(wdmain, "Data/")
wdassigndata <- paste0(wdproject, "Data/Assignment2groups/")
wdinputdata  <- paste0(wdproject, "Data/InputModel/")
wdoutputdata <- paste0(wdproject, "Data/OutputModel/")
wdfigs       <- paste0(wdproject, "Figures/")

### Ler arquivo d eentrada ----
input_file <- "inputdata_zmax0.1_Rlim2.5_Ma12.3_flag_good==1_MANGLE_logMstar_min10.5_dens_neighbors.csv"

gal <- fread(paste0(wddata, "environmental-quenching-data/InputModel/", input_file))
colnames(gal)[which(colnames(gal) == "AGN")] <- "AGN_char"
colnames(gal)[which(colnames(gal) == "SF_GSWLC")] <- "SF_char"

gal$SF_char[which(gal$SF_char == "Star forming")] <- "Star-forming"
gal$SF_char <- factor(gal$SF_char, levels = c("Star-forming", "Quiescent"))

gal$SF        <- ifelse(gal$SF_char == "Star-forming", 1, 0)
gal$AGN       <- ifelse(gal$AGN_char == "AGN", 1, 0)

data.model    <- subset(gal, gal$type == "Satellite")
data.centrais <- subset(gal, gal$type == "Central")

### Modelos ----

model1 <- glm(SF ~ 
                logMgroup +
                logMstar +
                logvelDisp_e +
                logRproj_rvir + 
                AGN,
              
              family = binomial(link = "logit"), 
              data = data.model)

model2 <- glm(SF ~ 
                poly(logMgroup,     1, raw = T) +
                poly(logMstar,      3, raw = T) +
                poly(logvelDisp_e,  3, raw = T) +
                poly(logRproj_rvir, 1, raw = T) + 
                AGN,
              
              family = binomial(link = "logit"), 
              data = data.model)

model3 <- glm(SF ~ 
                logdens_proj_Neq1 +
                logMstar +
                logvelDisp_e +
                AGN,
              
              family = binomial(link = "logit"), 
              data = data.model)

model4 <- glm(SF ~ 
                logdens_proj_Neq1 +
                poly(logMstar,      3, raw = T) +
                poly(logvelDisp_e,  3, raw = T) +
                AGN,
              
              family = binomial(link = "logit"), 
              data = data.model)

model5 <- glm(SF ~ 
                logdens_proj_Neq3 +
                poly(logMstar,      3, raw = T) +
                poly(logvelDisp_e,  3, raw = T) +
                AGN,
              
              family = binomial(link = "logit"), 
              data = data.model)

model6 <- glm(SF ~ 
                logdens_proj_Neq5 + 
                poly(logMstar,      3, raw = T) +
                poly(logvelDisp_e,  3, raw = T) +
                AGN,
              
              family = binomial(link = "logit"), 
              data = data.model)

model7 <- glm(SF ~ 
                logdens_n_Rproj0.5 +
                poly(logMstar,      3, raw = T) +
                poly(logvelDisp_e,  3, raw = T) +
                AGN,
              
              family = binomial(link = "logit"), 
              data = data.model)

model8 <- glm(SF ~ 
                logdens_n_Rproj1 + 
                poly(logMstar,      3, raw = T) +
                poly(logvelDisp_e,  3, raw = T) +
                AGN,
              
              family = binomial(link = "logit"), 
              data = data.model)

model9 <- glm(SF ~ 
                logdens_n_Rproj1 + 
                poly(logMstar,      3, raw = T) +
                poly(logvelDisp_e,  3, raw = T) +
                poly(logRproj_rvir, 1, raw = T) +
                AGN,
              
              family = binomial(link = "logit"), 
              data = data.model)

# Importância -----

imp1 <- as.data.frame(varImp(model1))
imp1 <- data.frame(overall = imp1$Overall, names = rownames(imp1))
imp1 <- imp1[order(imp1$overall, decreasing = T),]

imp2 <- as.data.frame(varImp(model2))
imp2 <- data.frame(overall = imp2$Overall, names = rownames(imp2))
imp2 <- imp2[order(imp2$overall, decreasing = T),]

imp3 <- as.data.frame(varImp(model3))
imp3 <- data.frame(overall = imp3$Overall, names = rownames(imp3))
imp3 <- imp3[order(imp3$overall, decreasing = T),]

imp4 <- as.data.frame(varImp(model4))
imp4 <- data.frame(overall = imp4$Overall, names = rownames(imp4))
imp4 <- imp4[order(imp4$overall, decreasing = T),]

imp5 <- as.data.frame(varImp(model5))
imp5 <- data.frame(overall = imp5$Overall, names = rownames(imp5))
imp5 <- imp5[order(imp5$overall, decreasing = T),]

imp6 <- as.data.frame(varImp(model6))
imp6 <- data.frame(overall = imp6$Overall, names = rownames(imp6))
imp6 <- imp6[order(imp6$overall, decreasing = T),]

imp7 <- as.data.frame(varImp(model7))
imp7 <- data.frame(overall = imp7$Overall, names = rownames(imp7))
imp7 <- imp7[order(imp7$overall, decreasing = T),]

imp8 <- as.data.frame(varImp(model8))
imp8 <- data.frame(overall = imp8$Overall, names = rownames(imp8))
imp8 <- imp8[order(imp8$overall, decreasing = T),]

imp9 <- as.data.frame(varImp(model9))
imp9 <- data.frame(overall = imp9$Overall, names = rownames(imp9))
imp9 <- imp9[order(imp9$overall, decreasing = T),]

# write.table(imp1, paste0("~/Work/Research/Astronomy/Projects/environmental-quenching/Data/Importance/", as.character(formula(model1))[3], ".txt"), row.names = F)
# write.table(imp2, paste0("~/Work/Research/Astronomy/Projects/environmental-quenching/Data/Importance/", as.character(formula(model2))[3], ".txt"), row.names = F)
# write.table(imp3, paste0("~/Work/Research/Astronomy/Projects/environmental-quenching/Data/Importance/", as.character(formula(model3))[3], ".txt"), row.names = F)
# write.table(imp4, paste0("~/Work/Research/Astronomy/Projects/environmental-quenching/Data/Importance/", as.character(formula(model4))[3], ".txt"), row.names = F)
# write.table(imp5, paste0("~/Work/Research/Astronomy/Projects/environmental-quenching/Data/Importance/", as.character(formula(model5))[3], ".txt"), row.names = F)
# write.table(imp6, paste0("~/Work/Research/Astronomy/Projects/environmental-quenching/Data/Importance/", as.character(formula(model6))[3], ".txt"), row.names = F)
# write.table(imp7, paste0("~/Work/Research/Astronomy/Projects/environmental-quenching/Data/Importance/", as.character(formula(model7))[3], ".txt"), row.names = F)
# write.table(imp8, paste0("~/Work/Research/Astronomy/Projects/environmental-quenching/Data/Importance/", as.character(formula(model8))[3], ".txt"), row.names = F)
# write.table(imp9, paste0("~/Work/Research/Astronomy/Projects/environmental-quenching/Data/Importance/", as.character(formula(model9))[3], ".txt"), row.names = F)

# Predictions for each model ----
predictions <- data.model
predictions$sf_b <- ifelse(predictions$SF == TRUE, 1, 0)

predictions$pred1 <- predict(model1, predictions, type = "response")
predictions$pred2 <- predict(model2, predictions, type = "response")
predictions$pred3 <- predict(model3, predictions, type = "response")
predictions$pred4 <- predict(model4, predictions, type = "response")
predictions$pred5 <- predict(model5, predictions, type = "response")
predictions$pred6 <- predict(model6, predictions, type = "response")
predictions$pred7 <- predict(model7, predictions, type = "response")
predictions$pred8 <- predict(model8, predictions, type = "response")
predictions$pred9 <- predict(model9, predictions, type = "response")

optimal1 <- optimalCutoff(predictions$sf_b, predictions$pred1)[1]
optimal2 <- optimalCutoff(predictions$sf_b, predictions$pred2)[1]
optimal3 <- optimalCutoff(predictions$sf_b, predictions$pred3)[1]
optimal4 <- optimalCutoff(predictions$sf_b, predictions$pred4)[1]
optimal5 <- optimalCutoff(predictions$sf_b, predictions$pred5)[1]
optimal6 <- optimalCutoff(predictions$sf_b, predictions$pred6)[1]
optimal7 <- optimalCutoff(predictions$sf_b, predictions$pred7)[1]
optimal8 <- optimalCutoff(predictions$sf_b, predictions$pred8)[1]
optimal9 <- optimalCutoff(predictions$sf_b, predictions$pred9)[1]

# Misclassifications ----

misClassError1 <- misClassError(predictions$sf_b, predictions$pred1, threshold = optimal1)
misClassError2 <- misClassError(predictions$sf_b, predictions$pred2, threshold = optimal2)
misClassError3 <- misClassError(predictions$sf_b, predictions$pred3, threshold = optimal3)
misClassError4 <- misClassError(predictions$sf_b, predictions$pred4, threshold = optimal4)
misClassError5 <- misClassError(predictions$sf_b, predictions$pred5, threshold = optimal5)
misClassError6 <- misClassError(predictions$sf_b, predictions$pred6, threshold = optimal6)
misClassError7 <- misClassError(predictions$sf_b, predictions$pred7, threshold = optimal7)
misClassError8 <- misClassError(predictions$sf_b, predictions$pred8, threshold = optimal8)
misClassError9 <- misClassError(predictions$sf_b, predictions$pred9, threshold = optimal9)

predictions$class1 <- ifelse(predictions$pred1 > optimal1, "Star forming", "Quiescent")
predictions$class2 <- ifelse(predictions$pred2 > optimal2, "Star-forming", "Quiescent")
predictions$class3 <- ifelse(predictions$pred3 > optimal3, "Star forming", "Quiescent")
predictions$class4 <- ifelse(predictions$pred4 > optimal4, "Star forming", "Quiescent")
predictions$class5 <- ifelse(predictions$pred5 > optimal5, "Star forming", "Quiescent")
predictions$class6 <- ifelse(predictions$pred6 > optimal6, "Star forming", "Quiescent")
predictions$class7 <- ifelse(predictions$pred7 > optimal7, "Star forming", "Quiescent")
predictions$class8 <- ifelse(predictions$pred8 > optimal8, "Star forming", "Quiescent")
predictions$class9 <- ifelse(predictions$pred9 > optimal9, "Star forming", "Quiescent")

# Model 1
predictions$misclass1 <- NA
predictions$misclass1[which(predictions$class1 == "Star forming" & predictions$SF_char == "Star forming")] <- "TP"
predictions$misclass1[which(predictions$class1 == "Quiescent"    & predictions$SF_char == "Quiescent")]    <- "TN"
predictions$misclass1[which(predictions$class1 == "Star forming" & predictions$SF_char == "Quiescent")]    <- "FP"
predictions$misclass1[which(predictions$class1 == "Quiescent"    & predictions$SF_char == "Star forming")] <- "FN"

# Model 2
predictions$misclass2 <- NA
predictions$misclass2[which(predictions$class2 == "Star forming" & predictions$SF_char == "Star forming")] <- "TP"
predictions$misclass2[which(predictions$class2 == "Quiescent"    & predictions$SF_char == "Quiescent")]    <- "TN"
predictions$misclass2[which(predictions$class2 == "Star forming" & predictions$SF_char == "Quiescent")]    <- "FP"
predictions$misclass2[which(predictions$class2 == "Quiescent"    & predictions$SF_char == "Star forming")] <- "FN"

# Model 3
predictions$misclass3 <- NA
predictions$misclass3[which(predictions$class3 == "Star forming" & predictions$SF_char == "Star forming")] <- "TP"
predictions$misclass3[which(predictions$class3 == "Quiescent"    & predictions$SF_char == "Quiescent")]    <- "TN"
predictions$misclass3[which(predictions$class3 == "Star forming" & predictions$SF_char == "Quiescent")]    <- "FP"
predictions$misclass3[which(predictions$class3 == "Quiescent"    & predictions$SF_char == "Star forming")] <- "FN"

# Model 4
predictions$misclass4 <- NA
predictions$misclass4[which(predictions$class4 == "Star forming" & predictions$SF_char == "Star forming")] <- "TP"
predictions$misclass4[which(predictions$class4 == "Quiescent"    & predictions$SF_char == "Quiescent")]    <- "TN"
predictions$misclass4[which(predictions$class4 == "Star forming" & predictions$SF_char == "Quiescent")]    <- "FP"
predictions$misclass4[which(predictions$class4 == "Quiescent"    & predictions$SF_char == "Star forming")] <- "FN"

# Model 5
predictions$misclass5 <- NA
predictions$misclass5[which(predictions$class5 == "Star forming" & predictions$SF_char == "Star forming")] <- "TP"
predictions$misclass5[which(predictions$class5 == "Quiescent"    & predictions$SF_char == "Quiescent")]    <- "TN"
predictions$misclass5[which(predictions$class5 == "Star forming" & predictions$SF_char == "Quiescent")]    <- "FP"
predictions$misclass5[which(predictions$class5 == "Quiescent"    & predictions$SF_char == "Star forming")] <- "FN"

# Model 6
predictions$misclass6 <- NA
predictions$misclass6[which(predictions$class6 == "Star forming" & predictions$SF_char == "Star forming")] <- "TP"
predictions$misclass6[which(predictions$class6 == "Quiescent"    & predictions$SF_char == "Quiescent")]    <- "TN"
predictions$misclass6[which(predictions$class6 == "Star forming" & predictions$SF_char == "Quiescent")]    <- "FP"
predictions$misclass6[which(predictions$class6 == "Quiescent"    & predictions$SF_char == "Star forming")] <- "FN"

# Model 7
predictions$misclass7 <- NA
predictions$misclass7[which(predictions$class7 == "Star forming" & predictions$SF_char == "Star forming")] <- "TP"
predictions$misclass7[which(predictions$class7 == "Quiescent"    & predictions$SF_char == "Quiescent")]    <- "TN"
predictions$misclass7[which(predictions$class7 == "Star forming" & predictions$SF_char == "Quiescent")]    <- "FP"
predictions$misclass7[which(predictions$class7 == "Quiescent"    & predictions$SF_char == "Star forming")] <- "FN"

# Model 8
predictions$misclass8 <- NA
predictions$misclass8[which(predictions$class8 == "Star forming" & predictions$SF_char == "Star forming")] <- "TP"
predictions$misclass8[which(predictions$class8 == "Quiescent"    & predictions$SF_char == "Quiescent")]    <- "TN"
predictions$misclass8[which(predictions$class8 == "Star forming" & predictions$SF_char == "Quiescent")]    <- "FP"
predictions$misclass8[which(predictions$class8 == "Quiescent"    & predictions$SF_char == "Star forming")] <- "FN"

# Model 9
predictions$misclass9 <- NA
predictions$misclass9[which(predictions$class9 == "Star forming" & predictions$SF_char == "Star forming")] <- "TP"
predictions$misclass9[which(predictions$class9 == "Quiescent"    & predictions$SF_char == "Quiescent")]    <- "TN"
predictions$misclass9[which(predictions$class9 == "Star forming" & predictions$SF_char == "Quiescent")]    <- "FP"
predictions$misclass9[which(predictions$class9 == "Quiescent"    & predictions$SF_char == "Star forming")] <- "FN"

write.csv(predictions, paste0(wdoutputdata, "predictions_allmodels.csv"), row.names = F)

# for (i in 1:nrow(FP)) {
#   # Images
#   path          <- "/home/muckler/Work/Research/Figures/FP_allmodels/"
#   file_name_img <- paste0(path, "specobjid_", FP$specobjid[i], "--RA_",FP$ra[i],"-DEC_",FP$dec[i],"_image.jpg")
#   file_name_spc <- paste0(path, "specobjid_", FP$specobjid[i], "--RA_",FP$ra[i],"-DEC_",FP$dec[i],"_spec.jpg")
#   
#   http_img      <- paste0("'http://skyserver.sdss.org/dr18/SkyserverWS/ImgCutout/getjpeg?ra=",FP$ra[i],"&dec=",FP$dec[i],"&scale=0.100000&width=512&height=512'")
#   http_spc      <- paste0("'http://skyserver.sdss.org/dr18/en/get/specById.ashx?ID=",FP$specobjid[i],"'")
#   
#   system(paste0("wget --output-document=",file_name_img, " ", http_img))
#   system(paste0("wget --output-document=",file_name_spc, " ", http_spc))
# }
# 
# for (i in 1:nrow(FN)) {
#   # Images
#   path          <- "/home/muckler/Work/Research/Figures/FN_allmodels/"
#   file_name_img <- paste0(path, "specobjid_", FN$specobjid[i], "--RA_",FN$ra[i],"-DEC_",FN$dec[i],"_image.jpg")
#   file_name_spc <- paste0(path, "specobjid_", FN$specobjid[i], "--RA_",FN$ra[i],"-DEC_",FN$dec[i],"_spec.jpg")
#   
#   http_img      <- paste0("'http://skyserver.sdss.org/dr18/SkyserverWS/ImgCutout/getjpeg?ra=",FN$ra[i],"&dec=",FN$dec[i],"&scale=0.100000&width=512&height=512'")
#   http_spc      <- paste0("'http://skyserver.sdss.org/dr18/en/get/specById.ashx?ID=",FN$specobjid[i],"'")
#   
#   system(paste0("wget --output-document=",file_name_img, " ", http_img))
#   system(paste0("wget --output-document=",file_name_spc, " ", http_spc))
# }