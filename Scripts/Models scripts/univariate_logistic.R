library(data.table)
library(InformationValue)
library(ggplot2)
library(scales)

## Lendo meus códigos ----
source("~/Work/Research/Astronomy/Projects/environmental-quenching/Scripts/Themes/my_theme.R")
source("~/Work/Research/Astronomy/Projects/environmental-quenching/Scripts/Themes/ggplot_theme_Publication-2.R")

zmax       <- 0.03
input_data <- paste0("inputdata_zmax",zmax,"_Rlim2.5_Ma12.3_flag_good==1_MANGLE_logMstar_min10.5.csv")

df    <- fread(paste0("~/Work/Research/Astronomy/Data/EnvQuenching/inputModel/GSWLC/", input_data))
df$SF <- ifelse(df$SF_GSWLC == "Star-forming", 1, 0)
df$SF <- as.factor(df$SF)

df    <- subset(df, df$type == "Satellite")

model_logMstar           <- glm(SF ~ logMstar, family = binomial(link = "logit"), data = df)
model_logvelDisp_e       <- glm(SF ~ logvelDisp_e, family = binomial(link = "logit"), data = df)
model_conc               <- glm(SF ~ conc, family = binomial(link = "logit"), data = df)
model_TType              <- glm(SF ~ TType, family = binomial(link = "logit"), data = df)

model_logRproj_rvir      <- glm(SF ~ logRproj_rvir, family = binomial(link = "logit"), data = df)
model_logMgroup          <- glm(SF ~ logMgroup, family = binomial(link = "logit"), data = df)
model_logdens_proj_Neq1  <- glm(SF ~ logdens_proj_Neq1, family = binomial(link = "logit"), data = df)
model_logdens_proj_Neq3  <- glm(SF ~ logdens_proj_Neq3, family = binomial(link = "logit"), data = df)
model_logdens_proj_Neq5  <- glm(SF ~ logdens_proj_Neq5, family = binomial(link = "logit"), data = df)
model_logdens_n_Rproj0.5 <- glm(SF ~ logdens_n_Rproj0.5, family = binomial(link = "logit"), data = df)
model_logdens_n_Rproj1   <- glm(SF ~ logdens_n_Rproj1, family = binomial(link = "logit"), data = df)

modelos <- list(model_logMstar,
                model_logvelDisp_e,
                model_conc,
                model_TType,
                model_logRproj_rvir,
                model_logMgroup,
                model_logdens_proj_Neq1,
                model_logdens_proj_Neq3,
                model_logdens_proj_Neq5,
                model_logdens_n_Rproj0.5,
                model_logdens_n_Rproj1)

for (i in 1:length(modelos)) {
  modelo <- modelos[[i]]
  prob_predita   <- predict(modelo, type = "response")
  optimal        <- optimalCutoff(df$SF, prob_predita)[1] 
  classe_predita <- ifelse(prob_predita > optimal, 1, 0)
  
  # Calcular a taxa de acerto
  acertos <- sum(classe_predita == df$SF)
  taxa_acerto <- acertos / nrow(df)
  taxa_acerto <- round(taxa_acerto*100, digits = 1)
  cat(paste0("Taxa de acerto do modelo ", i, ": ", taxa_acerto, "\n"))
}

# Figuras

df$pred<- predict(model_logMstar, type = "response")
ggplot(df, aes(x = logMstar, y = pred)) + 
  geom_line() +
  stat_smooth(formula = 'y ~ x', method = "glm", method.args = list(family = binomial)) + 
  scale_x_continuous(breaks = pretty_breaks(n = 10)) +
  scale_y_continuous(breaks = pretty_breaks(n = 10)) +
  xlab(mstar_label) + 
  ylab(fSFG_label) + 
  theme_Publication()

df$pred<- predict(model_TType, type = "response")
ggplot(df, aes(x = TType, y = pred)) + 
  geom_line() +
  stat_smooth(formula = 'y ~ x', method = "glm", method.args = list(family = binomial))

df$pred<- predict(model_logdens_proj_Neq5, type = "response")
ggplot(df, aes(x = logdens_proj_Neq5, y = pred)) + 
  geom_line() +
  stat_smooth(formula = 'y ~ x', method = "glm", method.args = list(family = binomial))

df$pred<- predict(model_logvelDisp_e, type = "response")
ggplot(df, aes(x = 10^logvelDisp_e, y = pred)) + 
  geom_line() +
  scale_x_log10(n.breaks = 15) + 
  stat_smooth(formula = 'y ~ x', method = "glm", method.args = list(family = binomial))

df$pred<- predict(model_logRproj_rvir, type = "response")
ggplot(df, aes(x = 10^logRproj_rvir, y = pred)) + 
  geom_line() +
  scale_x_log10(n.breaks = 15) + 
  scale_y_continuous(breaks = pretty_breaks(n = 10)) +
  stat_smooth(formula = 'y ~ x', method = "glm", method.args = list(family = binomial))

