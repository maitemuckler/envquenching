## Carregando pacotes ----
library(data.table)
library(dplyr)
library(ggplot2)
library(ggh4x)
library(scales)
library(latex2exp)

## Lendo meus códigos ----
source("~/Work/Research/Astronomy/Projects/environmental-quenching/Scripts/Themes/my_theme.R")
source("~/Work/Research/Astronomy/Projects/environmental-quenching/Scripts/Themes/ggplot_theme_Publication-2.R")

## Diretórios ----
wdmain       <- "~/Work/Research/Astronomy/Projects/environmental-quenching/"
wdcode       <- "~/Work/Research/Astronomy/Projects/environmental-quenching/Scripts/"
wdinputdata  <- "~/Work/Research/Astronomy/Data/EnvQuenching/inputModel/"
wdoutputdata <- "~/Work/Research/Astronomy/Data/EnvQuenching/outputModel/"
wdimportance <- "~/Work/Research/Astronomy/Data/EnvQuenching/importanceModel/"
wdfigs       <- "~/Work/Research/Astronomy/Projects/environmental-quenching/Figures/"

## Definindo opções ----
catalog <- "GSWLC"
critval <- 1.96 
zmax    <- 0.1
Rlim    <- 2.5
Ma      <- 12.3

tipo_de_modelo <- "NAGN"

if(zmax == 0.03){
  width_figs  <- 9
  height_figs <- 6
  size_text_facet <- 13
}

if(zmax == 0.1){
  width_figs  <- 9
  height_figs <- 7
  size_text_facet <- 12
}

## Definindo input e output files ----
bins_file  <- paste0("bins_", catalog, "_zmax", zmax, "_Ma", Ma, ".csv")
model_file <- paste0("model_", catalog, "_zmax", zmax, "_Ma", Ma, ".csv")
imp_file   <- paste0("importance_", catalog, "_zmax", zmax, "_Ma", Ma, ".csv")

## Lendo os dados ----
bins       <- fread(paste0(wdoutputdata, catalog, "/", tipo_de_modelo, "/", bins_file)) 
model      <- fread(paste0(wdoutputdata, catalog, "/", tipo_de_modelo, "/", model_file)) 
importance <- fread(paste0(wdimportance, catalog, "/", tipo_de_modelo, "/", imp_file))

p <- ggplot() + 
  geom_line(data = model, aes(x = 10^logRproj, y = pred, color = logMstar), linewidth = 1) + 
  geom_point(data = bins, aes(x = temp_rbins, y = fsf, color = logMstar), size = 2) +
  
  geom_errorbar(data = bins, aes(ymin = confint_lw, ymax = confint_up, x = temp_rbins, color = logMstar),
                width = 0.05, alpha = 0.4) + 
  
  geom_ribbon(data = model, aes(ymin = lwr, ymax = upr, x = 10^logRproj, fill = logMstar),
              inherit.aes = FALSE, alpha = 0.3) +
  
  # Centrais
  geom_point(data = bins, aes(x = 10^-1.15, y = fsf_c, color = logMstar), size = 3, shape = 1) +
  geom_errorbar(data = bins, aes(ymin = confint_lw_c, ymax = confint_up_c, x = 10^-1.15, color = logMstar),
                width = 0.05, alpha = 0.3) +
  
  facet_nested(logMgroup  ~ logvelDisp_e, labeller = label_parsed) + 
  
  scale_x_log10(n.breaks = 5) + 
  annotation_logticks(sides = "b") + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  coord_cartesian(xlim = 10^c(-1.28, 1.35), ylim = c(0,1)) + 
  
  scale_color_manual(values = c("#052d5a", "#d00000")) +
  scale_fill_manual(values = c("#052d5a", "#d00000")) +
  
  labs(x = expression(R/r["vir"]), 
       y = expression(f["SF"]), 
       color = expression(log['10']~M['★']),
       fill = expression(log['10']~M['★'])) + 
  
  theme_Publication() + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.position = "bottom",
        strip.text = element_text(size = size_text_facet),
        #strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) + 
  
  guides(fill = guide_legend(nrow = 1, byrow = F),
         color = guide_legend(nrow = 1, byrow = F),
         legend.spacing.y = unit(1.0, 'cm'))
p

ggsave(path = paste0(wdfigs, "Models/", catalog, "/", tipo_de_modelo, "/pdf/"),
       filename = paste0("logistica_Ma", Ma, "_zmax", zmax, "_", catalog, "_", tipo_de_modelo, ".pdf"),
       device = cairo_pdf, width = width_figs, height = height_figs, 
       units = "in", dpi = 600)

ggsave(path = paste0(wdfigs, "Models/", catalog, "/", tipo_de_modelo, "/png/"),
       filename = paste0("logistica_Ma", Ma, "_zmax", zmax, "_", catalog, "_", tipo_de_modelo, ".png"),
       width = width_figs, height = height_figs, 
       units = "in", dpi = 600)

# Importance ----

importance <- importance %>%
  arrange(overall)

if(tipo_de_modelo == "NAGN" & zmax == 0.1){
  expand <- 9
  superior_limit_importance <- 100
}

if(tipo_de_modelo == "NAGN" & zmax == 0.03){
  expand <- 3
  superior_limit_importance <- 25
}

if(tipo_de_modelo == "AGN" & zmax == 0.1){
  expand <- 0.6
}

if(tipo_de_modelo == "AGN" & zmax == 0.03){
  expand <- 0.2
}

mhalo_label <- TeX("$log_{10}M_{halo}$")
rvir_label  <- TeX("$log_{10}R/r_{vir}$")
mstar_label <- TeX("$log_{10}M_\\star$")
sigma_label <- TeX("$log_{10}\\sigma $")

importance$labels <- NA
importance$labels[which(importance$names == "logMstar")]      <- mstar_label
importance$labels[which(importance$names == "logMgroup")]     <- mhalo_label
importance$labels[which(importance$names == "logRproj_rvir")] <- rvir_label
importance$labels[which(importance$names == "logvelDisp_e")]  <- sigma_label

# Converter para data.frame
importance <- as.data.frame(importance)

# Reordenar o data frame
importance <- importance[order(-importance$overall), ]

# Atualizar os níveis do fator 'names' após a reordenação
importance$names <- factor(importance$names, levels = importance$names)

ggplot(importance, aes(x = reorder(names, -overall), y = overall)) +
  geom_col(fill = "skyblue3") + 
  geom_text(aes(label = round(overall, digits = 2)), size = 8, hjust = 0.5, vjust = -0.5) + 
  scale_y_continuous(breaks = pretty_breaks(), limits = c(0, superior_limit_importance)) + 
  scale_x_discrete(labels = importance$labels) +
  labs(y = "|z value|", x = "") + 
  #coord_flip() + 
  expand_limits(y = max(importance$overall) + expand) + 
  theme_Publication() + 
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 36),
        axis.text.x = element_text(angle = 20, hjust = 1),
        axis.text.y = element_text(margin = margin(r = 10)),
        legend.position = c(0.86, 0.88),
        legend.key.size = unit(1.5, "cm"),
        legend.spacing = unit(3, "cm"),
        legend.text = element_text(size = 24, margin = margin(l = 10)),
        legend.title = element_blank(),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.box.background = element_rect(colour = "black", linewidth = 2),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

ggsave(path = paste0(wdfigs, "Importance/", catalog, "/", tipo_de_modelo, "/pdf/"),
       filename = paste0("importance_Ma", Ma, "_zmax", zmax, "_", catalog, "_", tipo_de_modelo, ".pdf"),
       device = cairo_pdf, 
       width = 12, height = 7, 
       units = "in", dpi = 600)

ggsave(path = paste0(wdfigs, "Importance/", catalog, "/", tipo_de_modelo, "/png/"),
       filename = paste0("importance_Ma", Ma, "_zmax", zmax, "_", catalog, "_", tipo_de_modelo, ".png"),
       width = 12, height = 7, 
       units = "in", dpi = 600)


ggsave(path = paste0(wdfigs, "Importance/", catalog, "/", tipo_de_modelo, "/png/"),
       filename = paste0("importance_Ma", Ma, "_zmax", zmax, "_", catalog, "_", tipo_de_modelo, ".png"),
       width = 8, height = 11,
       units = "in", dpi = 600)
