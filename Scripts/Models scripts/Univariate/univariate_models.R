library(data.table)
library(InformationValue)
library(ggplot2)
library(scales)
library(dplyr)
library(binom)

## Lendo meus códigos ----
source("~/Work/Research/Astronomy/Projects/environmental-quenching/Scripts/Themes/my_theme.R")
source("~/Work/Research/Astronomy/Projects/environmental-quenching/Scripts/Themes/ggplot_theme_Publication-2.R")

figswd <- "~/Work/Research/Astronomy/Projects/environmental-quenching/Figures/Models/univariate_models"

critval         <- 1.96 
preditora       <- "logRproj_rvir"

if(preditora == "logvelDisp_e"){preditora_label <- label_logvelDisp_e; n_bins = 10; binagem = "cut"}
if(preditora == "logMstar"){preditora_label <- label_logMstar; n_bins = 10; binagem = "cut"}
if(preditora == "logRproj_rvir"){preditora_label <- label_logRproj_rvir; n_bins = 10; binagem = "cut"}
if(preditora == "TType"){preditora_label <- label_TType; n_bins = 10; binagem = "cut"}

input_data003 <- paste0("inputdata_zmax0.03_Rlim2.5_Ma12.3_flag_good==1_MANGLE_logMstar_min10.5.csv")
input_data01  <- paste0("inputdata_zmax0.1_Rlim2.5_Ma12.3_flag_good==1_MANGLE_logMstar_min10.5.csv")

colunas <- c("SF_GSWLC",
             "type",
             preditora)

df003 <- fread(paste0("~/Work/Research/Astronomy/Data/EnvQuenching/inputModel/GSWLC/", input_data003), select = colunas)
df01  <- fread(paste0("~/Work/Research/Astronomy/Data/EnvQuenching/inputModel/GSWLC/", input_data01), select = colunas)

colnames(df003)[which(colnames(df003) == preditora)] <- "preditora"
colnames(df01)[which(colnames(df01) == preditora)] <- "preditora"

df003 <- subset(df003, df003$type == "Satellite")
df01  <- subset(df01, df01$type == "Satellite")

df003$zmax <- 0.03
df01$zmax  <- 0.1

df003$SF <- ifelse(df003$SF_GSWLC == "Star-forming", 1, 0)
df003$SF <- as.factor(df003$SF)

df01$SF <- ifelse(df01$SF_GSWLC == "Star-forming", 1, 0)
df01$SF <- as.factor(df01$SF)

# Modelo
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

preddata003 <- with(df003, data.frame(preditora = seq(min(preditora), max(preditora), length = 100)))
preddata01  <- with(df01, data.frame(preditora = seq(min(preditora), max(preditora), length = 100)))

preddata003$accuracy <- taxa_acerto003
preddata01$accuracy  <- taxa_acerto01

preds003 <- predict(modelo003, newdata = preddata003, type = "link", se.fit = TRUE)
preds01  <- predict(modelo01, newdata = preddata01, type = "link", se.fit = TRUE)

# Intervalos de confiança

# zmax = 0.03
upr <- preds003$fit + (critval * preds003$se.fit)
lwr <- preds003$fit - (critval * preds003$se.fit)
fit <- preds003$fit

fit2 <- modelo003$family$linkinv(fit)
upr2 <- modelo003$family$linkinv(upr)
lwr2 <- modelo003$family$linkinv(lwr)

preddata003$lwr <- lwr2 
preddata003$upr <- upr2 

# zmax = 0.1
upr <- preds01$fit + (critval * preds01$se.fit)
lwr <- preds01$fit - (critval * preds01$se.fit)
fit <- preds01$fit

fit2 <- modelo01$family$linkinv(fit)
upr2 <- modelo01$family$linkinv(upr)
lwr2 <- modelo01$family$linkinv(lwr)

preddata01$lwr <- lwr2 
preddata01$upr <- upr2 

# labels
preddata003$zmax <- 0.03
preddata01$zmax  <- 0.1

# Dados binados

if(binagem == "cut"){
  
  bins_temp003 <- df003 %>%
    mutate(bins = cut(preditora, breaks = n_bins))
  
  bins_temp01 <- df01 %>%
    mutate(bins = cut(preditora, breaks = n_bins))
}

if(binagem == "ntile"){
  
  bins_temp003 <- df003 %>%
    mutate(bins = ntile(preditora, n_bins))
  
  bins_temp01 <- df01 %>%
    mutate(bins = ntile(preditora, n_bins))
}

bins_temp003 <- bins_temp003 %>%
  group_by(bins) %>%
  mutate(mediana_bins = median(preditora)) %>%
  summarize(
    mediana_bins = first(mediana_bins),
    total_linhas = n(),
    linhas_SF_igual_1 = sum(SF == 1),
    proporcao_SF_igual_1 = sum(SF == 1) / n()
  ) %>%
  ungroup()


bins_temp01 <- bins_temp01 %>%
  group_by(bins) %>%
  mutate(mediana_bins = median(preditora)) %>%
  summarize(
    mediana_bins = first(mediana_bins),
    total_linhas = n(),
    linhas_SF_igual_1 = sum(SF == 1),
    proporcao_SF_igual_1 = sum(SF == 1) / n()
  ) %>%
  ungroup()

bins_temp003$confint_up <- binom.confint(x = bins_temp003$linhas_SF_igual_1, 
                                         n = bins_temp003$total_linhas, 
                                         conf.level = 0.95, 
                                         method = "bayes", 
                                         type = "central")$upper
bins_temp003$confint_lw <- binom.confint(x = bins_temp003$linhas_SF_igual_1, 
                                         n = bins_temp003$total_linhas, 
                                         conf.level = 0.95, 
                                         method = "bayes", 
                                         type = "central")$lower

bins_temp01$confint_up <- binom.confint(x = bins_temp01$linhas_SF_igual_1, 
                                        n = bins_temp01$total_linhas, 
                                        conf.level = 0.95, 
                                        method = "bayes", 
                                        type = "central")$upper
bins_temp01$confint_lw <- binom.confint(x = bins_temp01$linhas_SF_igual_1, 
                                        n = bins_temp01$total_linhas, 
                                        conf.level = 0.95, 
                                        method = "bayes", 
                                        type = "central")$lower

bins_temp003$zmax <- 0.03
bins_temp01$zmax  <- 0.1

df        <- rbind(df003, df01)
bins_temp <- rbind(bins_temp003, bins_temp01)
preddata  <- rbind(preddata003, preddata01)

df$zmax        <- as.factor(df$zmax)
bins_temp$zmax <- as.factor(bins_temp$zmax)
preddata$zmax  <- as.factor(preddata$zmax)

rm(df003, df01, bins_temp003, bins_temp01, preddata003, preddata01)

g <- ggplot() + 
  geom_line(data = df, aes(x = preditora, y = pred, color = zmax)) +         
  geom_point(data = bins_temp, aes(y = proporcao_SF_igual_1, x = mediana_bins, color = zmax), size = 2) +
  geom_errorbar(data = bins_temp, aes(ymin = confint_lw, ymax = confint_up, x = mediana_bins, color = zmax), width = 0.025, alpha = 0.4) + 
  
  geom_ribbon(data = preddata, aes(ymin = lwr, ymax = upr, x = preditora, color = zmax, fill = zmax), inherit.aes = FALSE, alpha = 0.3) +
  #annotate("text", x = -Inf, y = 0, hjust = -0.1, vjust = 0, label = paste0("Accuracy = ", taxa_acerto003), color = "#8ac926", size = 5, family = "sans") + 
  #annotate("text", x = -Inf, y = 0, hjust = -0.1, vjust = -2, label = paste0("Accuracy = ", taxa_acerto01), color = "#ffb703", size = 5, family = "sans") + 
  
  geom_label(aes(x = -Inf, y = 0, label = paste0("Accuracy = ", preddata$accuracy[which(preddata$zmax == 0.03)][1], "%")), 
             fill = "white", color = "#8ac926", size = 5, family = "sans", label.size = 0, hjust = -0.1, vjust = 0) + 
  geom_label(aes(x = -Inf, y = 0, label = paste0("Accuracy = ", preddata$accuracy[which(preddata$zmax == 0.1)][1], "%")), 
             fill = "white", color = "#ffb703", size = 5, family = "sans", label.size = 0, hjust = -0.1, vjust = -1) + 
  
  scale_y_continuous(breaks = pretty_breaks(n = 5)) +
  scale_x_continuous(breaks = pretty_breaks(n = 5)) +
  scale_color_manual(values = c("#8ac926", "#ffb703")) + 
  scale_fill_manual(values = c("#8ac926", "#ffb703")) + 
  
  coord_cartesian(ylim = c(0, 1)) + 
  
  labs(x = preditora_label, y = fSFG_label, color = expression(z[max]), fill = expression(z[max])) + 
  
  theme_Publication() + 
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 34),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 34))

g

ggsave(path = figswd, plot = g,
       filename = paste0(preditora, ".pdf"),
       device = cairo_pdf, width = width_figs, height = height_figs, units = "in", dpi = 600)
