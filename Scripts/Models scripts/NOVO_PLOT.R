library(tidyverse)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggh4x)
library(scales)
library(binom)

critval <- 1.96 
`%!in%` = Negate(`%in%`)
source("Scripts/Themes/my_theme.R")

# Lendo os dados ----

gal <- fread("Data/input_model/input_model_GSWLC_zmax0.1_Ma12.3_flag_good==1_MANGLE.csv",
             select = c("logMgroup", 
                        "logMstar", 
                        "logvelDisp_e", 
                        "Rproj_rvir", 
                        "logRproj_rvir", 
                        "AGN", 
                        "type",
                        "SF_gswlc"))

colnames(gal)[which(colnames(gal) == "AGN")]      <- "AGN_char"
colnames(gal)[which(colnames(gal) == "SF_gswlc")] <- "SF_char"
gal$AGN_char <- as.factor(gal$AGN_char)
gal$SF_char  <- as.factor(gal$SF_char)

gal$SF          <- ifelse(gal$SF_char == "Star forming", 1, 0)
gal$AGN         <- ifelse(gal$AGN_char == "AGN", 1, 0)

data.satellites <- subset(gal, gal$type == "Satellite")
data.centrals   <- subset(gal, gal$type == "Central")

# Gerando modelo ----
fit_glm <- glm(SF ~ 
                 poly(logMgroup,     1, raw = T) +
                 poly(logMstar,      3, raw = T) +
                 poly(logvelDisp_e,  3, raw = T) +
                 poly(logRproj_rvir, 1, raw = T) + 
                 AGN,
               
               family = binomial(link = "logit"), 
               data = data.satellites)

# Fazendo a binagem ---- 

# Binagem - cuts
# data.satellites <- data.satellites %>% 
#   mutate(logMgroup_BIN  = cut_number(logMgroup,    n = 3, right = F)) %>%
#   mutate(logMstar_BIN   = cut_number(logMstar,     n = 2, right = F)) %>%
#   mutate(logvelDisp_BIN = cut_number(logvelDisp_e, n = 3, right = F)) %>%
#   mutate(Rproj_rvir_BIN = cut_number(Rproj_rvir,   n = 5, right = F))

data.satellites <- data.satellites %>% 
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

#c(0, quantile(data.satellites$Rproj_rvir, seq(0, 1, by = 0.2)))

# Número de linhas por quadrante
data.satellites_countquad <- data.satellites %>% 
  count(logMgroup_BIN, logMstar_BIN, logvelDisp_BIN, AGN_char)

# Número de linhas por bin
data.satellites_countbin <- data.satellites %>% 
  count(logMgroup_BIN, logMstar_BIN, logvelDisp_BIN, AGN_char, Rproj_rvir_BIN)

# Número de linhas por classe (SF ou Q) em cada bin
data.satellites_countbinClass <- data.satellites %>% 
  count(logMgroup_BIN, logMstar_BIN, logvelDisp_BIN, AGN_char, Rproj_rvir_BIN, SF_char)

# Para calcular todas as medianas e fsf, vamos pegar uma linha de cada tabela data.satellites_countbin

tabelona <- data.frame(Rproj_rvir_BIN = character(0),
                       median_logMgroup = numeric(0),
                       median_logMstar = numeric(0),
                       median_logvelDisp = numeric(0),
                       median_logRproj_rvir = numeric(0),
                       AGN = numeric(0),
                       logMgroup_BIN = character(0),
                       logMstar_BIN = character(0),   
                       logvelDisp_BIN = character(0),
                       AGN_char = character(0),
                       n.x = numeric(0),
                       n.y = numeric(0),                 
                       fsf = numeric(0),
                       pred = numeric(0),
                       lwr = numeric(0),
                       upr = numeric(0),             
                       confint_up = numeric(0),
                       confint_lw = numeric(0))
i=1
j=1
k=1
l=1

for (i in 1:3) {
  for (j in 1:2) {
    for (k in 1:3) {
      for (l in 1:2) {
        
        # QUADRANTE 1
        bin_logMgroup  <- levels(data.satellites_countquad$logMgroup_BIN)[i]
        bin_logMstar   <- levels(data.satellites_countquad$logMstar_BIN)[j]
        bin_logvelDisp <- levels(data.satellites_countquad$logvelDisp_BIN)[k]
        bin_AGN        <- levels(data.satellites_countquad$AGN_char)[l]
        
        # Dados do quadrante 1
        data.satellites_quad <- data.satellites %>%
          filter(logMgroup_BIN  == bin_logMgroup &
                   logMstar_BIN   == bin_logMstar &
                   logvelDisp_BIN == bin_logvelDisp &
                   AGN_char       == bin_AGN)
        
        # Valores de mediana para o modelo do quadrante 1
        # para cada bin de Rproj, qual a mediana nas outras variáveis)
        
        data.satellites_quad_model <- data.satellites_quad %>%
          group_by(Rproj_rvir_BIN) %>% 
          summarise(median_logMgroup = median(logMgroup), 
                    median_logMstar = median(logMstar),
                    median_logvelDisp = median(logvelDisp_e),
                    median_logRproj_rvir = median(logRproj_rvir))
        
        # Adicionando coluna de AGN que estou usando no momento:
        data.satellites_quad_model$AGN <- ifelse(bin_AGN == "AGN", 1, 0)
        
        preds <- predict(fit_glm, 
                         newdata = data.satellites_quad_model %>% 
                           select(-Rproj_rvir_BIN) %>%
                           rename("logMgroup" = "median_logMgroup",
                                  "logMstar" = "median_logMstar",
                                  "logvelDisp_e" = "median_logvelDisp",
                                  "logRproj_rvir" = "median_logRproj_rvir"),
                         type = "link", se.fit = TRUE)
        
        pred <- fit_glm$family$linkinv(preds$fit)
        upr  <- fit_glm$family$linkinv(preds$fit + (critval * preds$se.fit))
        lwr  <- fit_glm$family$linkinv(preds$fit - (critval * preds$se.fit))
        
        
        # Pontos dos dados:
        index <- which(data.satellites_countbin$logMgroup_BIN == bin_logMgroup &
                         data.satellites_countbin$logMstar_BIN == bin_logMstar &
                         data.satellites_countbin$logvelDisp_BIN == bin_logvelDisp &
                         data.satellites_countbin$AGN_char == bin_AGN)
        
        data.satellites_quad_bin <- data.satellites_countbin[index,]
        dados_plot <- merge(data.satellites_quad_model, data.satellites_quad_bin, by = "Rproj_rvir_BIN")
        
        data.satellites_countbinSF <- data.satellites_countbinClass %>%
          filter(SF_char == "Star forming")
        
        aux <- data.satellites_countbinSF %>% 
          filter(logMgroup_BIN == bin_logMgroup &
                   logMstar_BIN == bin_logMstar &
                   logvelDisp_BIN == bin_logvelDisp &
                   AGN_char == bin_AGN)
        
        # tem dados para todos os bins de Rproj_rvir?
        if(nrow(aux) == length(levels(data.satellites$Rproj_rvir_BIN))){
          
        }else{
          falta <- levels(data.satellites$Rproj_rvir_BIN)[which(
            levels(data.satellites$Rproj_rvir_BIN) %!in% unique(aux$Rproj_rvir_BIN))]
          
          aux <- aux %>%
            add_row(logMgroup_BIN = bin_logMgroup,
                    logMstar_BIN = bin_logMstar,
                    logvelDisp_BIN = bin_logvelDisp,
                    AGN_char = bin_AGN,
                    Rproj_rvir_BIN = falta,
                    SF_char = "Star forming",
                    n = 0)
        }
        
        dados_plot <- merge(dados_plot, aux %>% select(Rproj_rvir_BIN, n), 
                            by = "Rproj_rvir_BIN") # cuidado aqui
        
        dados_plot <- dados_plot %>%
          mutate(fsf = n.y/n.x)
        
        dados_plot$pred <- pred
        dados_plot$lwr  <- lwr
        dados_plot$upr  <- upr
        
        dados_plot$confint_up <- binom.confint(x = dados_plot$n.y, n = dados_plot$n.x, 
                                               conf.level = 0.95, 
                                               method = "bayes", 
                                               type = "central")$upper
        dados_plot$confint_lw <- binom.confint(x = dados_plot$n.y, n = dados_plot$n.x, 
                                               conf.level = 0.95, 
                                               method = "bayes", 
                                               type = "central")$lower
        
        tabelona <- rbind(tabelona, dados_plot)
        
      }
    }
  }
}

tabelona$label <- paste0(tabelona$logMstar_BIN,'-',tabelona$AGN_char)

width_figs  <- 9
height_figs <- 7
size_text_facet <- 12

ggplot(tabelona, aes(x = 10^(median_logRproj_rvir), color = label)) + 
  geom_path(aes(y = pred, color = label, linetype = AGN_char), linewidth = 1) + 
  geom_point(aes(y = fsf, color = label, shape = AGN_char), size = 2) + 
  
  geom_errorbar(aes(ymin = confint_lw, ymax = confint_up, x = 10^(median_logRproj_rvir), 
                    color = label),width = 0.05, alpha = 0.4) + 
  
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = 10^(median_logRproj_rvir), 
                  fill = label), inherit.aes = FALSE, alpha = 0.3) +
  
  facet_nested(logMgroup_BIN  ~ logvelDisp_BIN) + 
  scale_x_log10(n.breaks = 5) + 
  annotation_logticks(sides = "b") + 
  #scale_y_continuous(breaks = pretty_breaks(n = 5)) + 
  coord_cartesian(xlim = 10^c(-1.28, 1.35), ylim = c(0,1)) + 
  labs(x = expression(R/r["vir"]), 
       y = expression(f["SF"])) + 
  
  scale_color_manual(values = c("#052d5a", "#052d5a", "#d00000", "#d00000")) +
  scale_fill_manual(values = c("#052d5a", "#052d5a", "#d00000", "#d00000")) +
  scale_linetype_manual(values = c("solid","dotdash")) + 
  scale_shape_manual(values = c(19, 15)) + 
  
  labs(x = expression(R/r["vir"]), 
       y = expression(f["SF"]), 
       color = expression(log['10']~M['★']),
       fill = expression(log['10']~M['★']),
       shape = "AGN",
       linetype = "AGN") + 
  my_theme + 
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

ggsave(filename = "medianadosbins.pdf",
       device = cairo_pdf, width = width_figs, height = height_figs, 
       units = "in", dpi = 600)

ggsave(filename = "medianadosbins.png",
       width = width_figs, height = height_figs, 
       units = "in", dpi = 600)



#ggsave(filename = paste0(titulo, ".png"))
