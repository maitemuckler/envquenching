library(caret)
library(tidyr)
library(ggplot2)
library(scales)

set.seed(123)

source("~/Work/Research/Astronomy/Projects/environmental-quenching/Scripts/Themes/my_theme.R")

load("~/Work/Research/Astronomy/Projects/environmental-quenching/Data/satelites.RData")

data.satellites$SF               <- as.factor(data.satellites$SF)
data.satellites$AGN              <- as.factor(data.satellites$AGN)
data.satellites$SF_GSWLC_central <- as.factor(data.satellites$SF_GSWLC_central)
data.satellites$AGN_central      <- as.factor(data.satellites$AGN_central)

trainIndex <- createDataPartition(data.satellites$SF, p = 0.7, list = FALSE)
train_data <- data.satellites[trainIndex, ]
test_data  <- data.satellites[-trainIndex, ]

rm(data.satellites)
rm(trainIndex)

# Ajustando o modelo de regressão logística (internas) ----------------------------------------------

variaveis <- c("logvelDisp_e",
               "logMstar",
               "vlos_vvir",
               "B_T_r",
               "TType",
               "conc",
               "d4000_n",
               "e")

facet_labels <- labeller(
  logvelDisp_e = sigmagal_label,
  logMstar = mstar_label,
  vlos_vvir = vlosvvir,
  B_T_r = B/T,
  TType = T-Type,
  conc = Concentration,
  d4000_n = Dn4000,
  e = Ellipticity)

# Modelos simples sem variáveis ambiental:
formula1 <- as.formula(paste("SF ~ logvelDisp_e"))
formula2 <- as.formula(paste("SF ~ logMstar"))
formula3 <- as.formula(paste("SF ~ vlos_vvir"))
formula4 <- as.formula(paste("SF ~ B_T_r"))
formula5 <- as.formula(paste("SF ~ TType"))
formula6 <- as.formula(paste("SF ~ conc"))
formula7 <- as.formula(paste("SF ~ d4000_n"))
formula8 <- as.formula(paste("SF ~ e"))

modelo1 <- glm(formula1, family = binomial(link = "logit"), data = train_data)
modelo2 <- glm(formula2, family = binomial(link = "logit"), data = train_data)
modelo3 <- glm(formula3, family = binomial(link = "logit"), data = train_data)
modelo4 <- glm(formula4, family = binomial(link = "logit"), data = train_data)
modelo5 <- glm(formula5, family = binomial(link = "logit"), data = train_data)
modelo6 <- glm(formula6, family = binomial(link = "logit"), data = train_data)
modelo7 <- glm(formula7, family = binomial(link = "logit"), data = train_data)
modelo8 <- glm(formula8, family = binomial(link = "logit"), data = train_data)

test_data$pred_prob1 <- predict(modelo1, newdata = test_data, type = "response")
test_data$pred_prob2 <- predict(modelo2, newdata = test_data, type = "response")
test_data$pred_prob3 <- predict(modelo3, newdata = test_data, type = "response")
test_data$pred_prob4 <- predict(modelo4, newdata = test_data, type = "response")
test_data$pred_prob5 <- predict(modelo5, newdata = test_data, type = "response")
test_data$pred_prob6 <- predict(modelo6, newdata = test_data, type = "response")
test_data$pred_prob7 <- predict(modelo7, newdata = test_data, type = "response")
test_data$pred_prob8 <- predict(modelo8, newdata = test_data, type = "response")

optimal1 <- optimalCutoff(test_data$SF, test_data$pred_prob1)
optimal2 <- optimalCutoff(test_data$SF, test_data$pred_prob2)
optimal3 <- optimalCutoff(test_data$SF, test_data$pred_prob3)
optimal4 <- optimalCutoff(test_data$SF, test_data$pred_prob4)
optimal5 <- optimalCutoff(test_data$SF, test_data$pred_prob5)
optimal6 <- optimalCutoff(test_data$SF, test_data$pred_prob6)
optimal7 <- optimalCutoff(test_data$SF, test_data$pred_prob7)
optimal8 <- optimalCutoff(test_data$SF, test_data$pred_prob8)

test_data$classe_predita1 <- ifelse(test_data$pred_prob1 > optimal1, 1, 0)
test_data$classe_predita2 <- ifelse(test_data$pred_prob2 > optimal2, 1, 0)
test_data$classe_predita3 <- ifelse(test_data$pred_prob3 > optimal3, 1, 0)
test_data$classe_predita4 <- ifelse(test_data$pred_prob4 > optimal4, 1, 0)
test_data$classe_predita5 <- ifelse(test_data$pred_prob5 > optimal5, 1, 0)
test_data$classe_predita6 <- ifelse(test_data$pred_prob6 > optimal6, 1, 0)
test_data$classe_predita7 <- ifelse(test_data$pred_prob7 > optimal7, 1, 0)
test_data$classe_predita8 <- ifelse(test_data$pred_prob8 > optimal8, 1, 0)

# Calcular a acurácia para cada modelo
accuracy1 <- 1 - misClassError(test_data$SF, test_data$pred_prob1, threshold = optimal1)
accuracy2 <- 1 - misClassError(test_data$SF, test_data$pred_prob2, threshold = optimal2)
accuracy3 <- 1 - misClassError(test_data$SF, test_data$pred_prob3, threshold = optimal3)
accuracy4 <- 1 - misClassError(test_data$SF, test_data$pred_prob4, threshold = optimal4)
accuracy5 <- 1 - misClassError(test_data$SF, test_data$pred_prob5, threshold = optimal5)
accuracy6 <- 1 - misClassError(test_data$SF, test_data$pred_prob6, threshold = optimal6)
accuracy7 <- 1 - misClassError(test_data$SF, test_data$pred_prob7, threshold = optimal7)
accuracy8 <- 1 - misClassError(test_data$SF, test_data$pred_prob8, threshold = optimal8)

# plot_data <- data.frame(
#   logvelDisp_e = train_data$logvelDisp_e,
#   logMstar = train_data$logMstar,
#   vlos_vvir = train_data$vlos_vvir,
#   B_T_r = train_data$B_T_r,
#   TType = train_data$TType,
#   conc = train_data$conc,
#   d4000_n = train_data$d4000_n,
#   e = train_data$e,
#   SF_char = train_data$SF_char
# )

variaveis <- c("logvelDisp_e",
               "logMstar",
               "vlos_vvir",
               "B_T_r",
               "TType",
               "conc")

plot_data <- data.frame(
  logvelDisp_e = train_data$logvelDisp_e,
  logMstar = train_data$logMstar,
  vlos_vvir = train_data$vlos_vvir,
  B_T_r = train_data$B_T_r,
  TType = train_data$TType,
  conc = train_data$conc,
  SF_char = train_data$SF_char
)

plot_data <- plot_data[sample(1:nrow(plot_data), size = 0.25*nrow(plot_data)),]

plot_data_long <- pivot_longer(plot_data, cols = logvelDisp_e:conc, names_to = "Variable", values_to = "Value")

accuracy_data <- data.frame(
  Variable = variaveis,
  Accuracy = c(accuracy1, accuracy2, accuracy3, accuracy4, accuracy5, accuracy6)
)

# plot_data_long$Accuracy <- NA
# for (nm_variavel in variaveis) {
#   acuracia <- accuracy_data$Accuracy[which(accuracy_data$Variable == nm_variavel)]
#   plot_data_long$Accuracy[which(plot_data_long$Variable == nm_variavel)] <- acuracia
# }

plot_data_long$Variable <- factor(plot_data_long$Variable, 
                                  levels = accuracy_data$Variable[order(-accuracy_data$Accuracy)])

ggplot(plot_data_long, aes(x = Value, y = as.numeric(SF_char == "Star-forming"))) +
  geom_point(alpha = 0.1, color = "gray") +  # Pontos dos dados observados
  stat_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE, color = "black", linewidth = 2) +  
  scale_x_continuous(breaks = pretty_breaks()) + 
  facet_wrap(. ~ factor(Variable), scales = "free_x", ncol = 3) + 
  labs(x = "", y = fSFG_label) +
  geom_text(data = accuracy_data, aes(x = Inf, y = 0.95, 
                                      label = paste0("Accuracy: ", round(Accuracy, 2))),
            hjust = 1.1, vjust = 1.5, color = "black", size = 5) + 
  theme_Publication() + 
  theme(strip.background = element_rect(colour = "black", fill = "#caebc8"),
        strip.text = element_text(size = 12, face = "plain"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))

ggsave(path = wdfigs,
       filename = paste0("models_1.pdf"),
       device = cairo_pdf, width = 16, height = 8, 
       units = "in", dpi = 600)

ggsave(path = paste0(wdfigs, "png/"),
       filename = paste0("models_1.png"),
       width = 16, height = 8, 
       units = "in", dpi = 600)

# Ajustando o modelo de regressão logística (externas) ----------------------------------------------
variaveis <- c("logMgroup",
               "logRproj_rvir",
               "logdens_n_Rproj0.5",
               "logdens_n_Rproj1",
               "logdens_proj_Neq1",
               "logdens_proj_Neq3",
               "logdens_proj_Neq5")

modelos <- data.frame(Variable = variaveis)
modelos$Accuracy <- NA

for (i in 1:length(variaveis)) {
  variavel            <- modelos$Variable[i]
  formula             <- as.formula(paste("SF ~", variavel))
  modelo              <- glm(formula, family = binomial(link = "logit"), data = train_data)
  pred_prob           <- predict(modelo, newdata = test_data, type = "response")
  optimal             <- optimalCutoff(test_data$SF, pred_prob)
  classe_predita      <- ifelse(pred_prob > optimal, 1, 0)
  modelos$Accuracy[i] <- 1 - misClassError(test_data$SF, pred_prob, threshold = optimal)
}

plot_data <- train_data %>%
  select(c(variaveis, "SF_char"))

plot_data      <- plot_data[sample(1:nrow(plot_data), size = 0.25*nrow(plot_data)),]
plot_data_long <- pivot_longer(plot_data, cols = variaveis[1]:variaveis[length(variaveis)], names_to = "Variable", values_to = "Value")

plot_data_long$Variable <- factor(plot_data_long$Variable, 
                                  levels = modelos$Variable[order(-modelos$Accuracy)])

ggplot(plot_data_long, aes(x = Value, y = as.numeric(SF_char == "Star-forming"))) +
  geom_point(alpha = 0.1, color = "gray") + 
  stat_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE, color = "black", linewidth = 2) +  
  scale_x_continuous(breaks = pretty_breaks()) + 
  facet_wrap(. ~ factor(Variable), scales = "free_x", ncol = 3) + 
  labs(x = "", y = fSFG_label) +
  geom_text(data = modelos, aes(x = Inf, y = 0.95, 
                                      label = paste0("Accuracy: ", round(Accuracy, 2))),
            hjust = 1.1, vjust = 1.5, color = "black", size = 5) + 
  theme_Publication() + 
  theme(strip.background = element_rect(colour = "black", fill = "#caebc8"),
        strip.text = element_text(size = 12, face = "plain"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))

ggsave(path = wdfigs,
       filename = paste0("models_2.pdf"),
       device = cairo_pdf, width = 16, height = 8, 
       units = "in", dpi = 600)

ggsave(path = paste0(wdfigs, "png/"),
       filename = paste0("models_2.png"),
       width = 16, height = 8, 
       units = "in", dpi = 600)
