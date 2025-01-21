library(lme4)
library(InformationValue)

dados <- data.frame(SF = data.model$SF,
                    groupID = data.model$groupID,
                    AGN = data.model$AGN,
                    logMstar = scale(data.model$logMstar, center = TRUE, scale = TRUE),
                    logMgroup = scale(data.model$logMgroup, center = TRUE, scale = TRUE),
                    logvelDisp_e = scale(data.model$logvelDisp_e, center = TRUE, scale = TRUE))

dados$groupID <- as.factor(dados$groupID)
#dados$SF <- as.factor(dados$SF)
dados$AGN <- as.factor(dados$AGN)

m1 <- glmer(SF ~ logMstar + logMgroup + logvelDisp_e + AGN + (1 | groupID),
            data = dados,
            family = binomial)
summary(m1)

test <- dados
test$sf_b <- ifelse(test$SF == TRUE, 1, 0)
test$pred1 <- predict(m1, test, type = "response")
optimal1 <- optimalCutoff(test$sf_b, test$pred1)[1]
misClassError(test$sf_b, test$pred1, threshold = optimal1) #0.2114

# Dados desbalanceados

barplot(prop.table(table(dados$AGN)),
        col = rainbow(2),
        ylim = c(0, 1),
        main = "Class Distribution")

library(caret)
set.seed(1234)

index <- createDataPartition(dados$AGN, p=0.8, list=FALSE)
train <- dados[index,]
test  <- dados[-index,]

which(colnames(dados) == "AGN")

trainup <- upSample(x=train[,-3], y=train$AGN)

table(trainup$Class)

colnames(trainup)[which(colnames(trainup) == "Class")] <- "AGN"

m2 <- glmer(SF ~ logMstar + logMgroup + logvelDisp_e + AGN + (1 | groupID),
            data = trainup,
            family = binomial)
summary(m2)

test <- trainup
test$sf_b <- ifelse(test$SF == TRUE, 1, 0)
test$pred1 <- predict(m2, test, type = "response")
optimal1 <- optimalCutoff(test$sf_b, test$pred1)[1]
misClassError(test$sf_b, test$pred1, threshold = optimal1) #0.0755

traindown <- downSample(x=train[,-3], y=train$AGN, yname = "AGN")
m3 <- glmer(SF ~ logMstar + logMgroup + logvelDisp_e + AGN + (1 | groupID),
            data = traindown,
            family = binomial)

test <- traindown
test$sf_b <- ifelse(test$SF == TRUE, 1, 0)
test$pred1 <- predict(m3, test, type = "response")
optimal1 <- optimalCutoff(test$sf_b, test$pred1)[1]
misClassError(test$sf_b, test$pred1, threshold = optimal1) #0.18

# 

m4 <- glmer(SF ~ logMstar + logMgroup + logvelDisp_e + (1 | groupID),
            data = dados,
            family = binomial)

test <- dados
test$sf_b <- ifelse(test$SF == TRUE, 1, 0)
test$pred1 <- predict(m4, test, type = "response")
optimal1 <- optimalCutoff(test$sf_b, test$pred1)[1]
misClassError(test$sf_b, test$pred1, threshold = optimal1) #0.0755
