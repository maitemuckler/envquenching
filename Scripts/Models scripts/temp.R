### Bibliotecas ----
library(ggplot2)
library(scales)
library(recipes)
library(data.table)
library(viridis)

### Definir diretórios ----
wdmain       <- "~/Work/Research/Astronomy/"
wdproject    <- paste0(wdmain, "Projects/environmental-quenching/")
wdcode       <- paste0(wdproject, "Scripts/")
wddata       <- paste0(wdmain, "Data/")
wdassigndata <- paste0(wdproject, "Data/Assignment2groups/")
wdinputdata  <- paste0(wdproject, "Data/InputModel/")
wdoutputdata <- paste0(wdproject, "Data/OutputModel/")
wdfigs       <- paste0(wdproject, "Figures/")

### Ler os dados ----
predictions <- fread(paste0(wdoutputdata, "predictions_allmodels.csv"))
predictions <- fread("~/Work/Research/Astronomy/Data/environmental-quenching-data/")
### Discretizando variáveis ----

breaks <- 30

breaks_logMstar           <- pretty(range(predictions$logMstar), n = breaks) # Obter os limites dos intervalos
midpoints_logMstar        <- (head(breaks_logMstar, -1) + tail(breaks_logMstar, -1)) / 2 # Calcular os pontos médios dos intervalos
predictions$logMstar_cut  <- cut(predictions$logMstar, breaks = breaks_logMstar, labels = midpoints_logMstar)

breaks_logvelDisp          <- pretty(range(predictions$logvelDisp_e), n = breaks)
midpoints_logvelDisp       <- (head(breaks_logvelDisp, -1) + tail(breaks_logvelDisp, -1)) / 2 
predictions$logvelDisp_cut <- cut(predictions$logvelDisp_e, breaks = breaks_logvelDisp, labels = midpoints_logvelDisp)

breaks_logsSFR_GSWLC          <- pretty(range(predictions$logsSFR_GSWLC), n = breaks)
midpoints_logsSFR_GSWLC       <- (head(breaks_logsSFR_GSWLC, -1) + tail(breaks_logsSFR_GSWLC, -1)) / 2 
predictions$logsSFR_GSWLC_cut <- cut(predictions$logsSFR_GSWLC, breaks = breaks_logsSFR_GSWLC, labels = midpoints_logsSFR_GSWLC)

predictions$P_edge_on_cut <- ifelse(predictions$P_edge_on > 0.9, "P_edge_on > 0.9", "P_edge_on <= 0.9")
predictions$e_cut        <- cut(predictions$e, breaks = 4)

### Criando factors ----
predictions$SF_char[which(predictions$SF_char == "Star forming")] <- "Star-forming"
predictions$SF_char <- factor(predictions$SF_char, levels = c("Star-forming", "Quiescent"))

### Tabela com misclassified ----

# Mesmo resultado em todos os modelos:
same_class <- predictions %>% filter(if_all(misclass2:misclass9, ~ . == misclass1))
same_class <- same_class %>%
  select(-paste0("misclass", 2:9))
colnames(same_class)[which(colnames(same_class) == "misclass1")] <- "misclass"

corretas <- same_class[which(same_class$misclass == "TP" | same_class$misclass == "TN"),]
erradas  <- same_class[which(same_class$misclass == "FP" | same_class$misclass == "FN"),]

# Gráficos ----


df_means <- aggregate(e ~ logMstar_cut + logvelDisp_cut + AGN_char, data = same_class, mean)

ggplot(df_means, aes(x = logMstar_cut, y = logvelDisp_cut, fill = e)) +
  geom_tile() + 
  scale_fill_viridis(option = "inferno", direction = -1) + 
  facet_wrap( . ~ AGN_char) + 
  scale_x_discrete(breaks = midpoints_logMstar[seq(1, length(midpoints_logMstar), by = 4)], 
                   labels = round(midpoints_logMstar[seq(1, length(midpoints_logMstar), by = 4)], digits = 2)) + 
  scale_y_discrete(breaks = midpoints_logvelDisp[seq(1, length(midpoints_logvelDisp), by = 4)],
                   labels = round(midpoints_logvelDisp[seq(1, length(midpoints_logvelDisp), by = 4)], digits = 2)) + 
  theme_bw()


# -------------

ggplot() + 
  geom_point(data = erradas, aes(x = logMstar, y = logsSFR_GSWLC, color = misclass), alpha = 0.5) + 
  geom_density_2d(data = predictions, aes(x = logMstar, y = logsSFR_GSWLC, color = SF_char), linewidth = 1) + 
  scale_color_manual(values = values) + 
  #scale_alpha_manual(values = alphas) + 
  scale_y_continuous(breaks = pretty_breaks(n=10)) + 
  scale_x_continuous(breaks = pretty_breaks(n=10)) + 
  labs(x = expression(log['10']~M['★']),
       y = expression(log['10']~sSFR), 
       color = "Classification",
       alpha = "Classification",
       size = "Classification") + 
  theme_bw() + 
  my_theme + 
  my_guides


ggplot() + 
  geom_point(data = erradas, aes(x = logMstar, y = logvelDisp_e, color = misclass), alpha = 0.5) + 
  geom_density_2d(data = predictions, aes(x = logMstar, y = logvelDisp_e, color = misclass), size = 1) + 
  scale_color_manual(values = cores) + 
  #scale_alpha_manual(values = alphas) + 
  scale_y_continuous(breaks = pretty_breaks(n=10)) + 
  scale_x_continuous(breaks = pretty_breaks(n=10)) + 
  labs(x = expression(log['10']~M['★']),
       y = expression(log['10']~sigma), 
       color = "Classification",
       alpha = "Classification",
       size = "Classification") + 
  theme_bw() + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.position = "bottom",
        strip.text = element_text(size = 13)) + 
  guides(fill = guide_legend(nrow = 2, byrow = F),
         color = guide_legend(nrow = 2, byrow = F),
         linetype = guide_legend(nrow = 2, byrow = F),
         shape = guide_legend(nrow = 2, byrow = F),
         legend.spacing.y = unit(1.0, 'cm'))



ggplot(predictions) + 
  geom_point(aes(x = logMstar, y = logvelDisp_e, color = e), alpha = 0.2)

# elipticidade do bojo (Bulge ellipticity)
ggplot(same_class) + 
  geom_density(aes(x = e, color = misclass))

ggplot(same_class) + 
  geom_density(aes(x = e, color = misclass), linewidth = op_linewidth) + 
  facet_grid(P_edge_on_cut ~ .) + 
  scale_x_continuous(breaks = pretty_breaks(n = 10)) + 
  scale_colour_manual(values = cores) + 
  theme_bw() + 
  theme(legend.position = "bottom") 


# nb
g2 <- ggplot(same_class) + 
  geom_density(aes(x = nb, color = misclass1), linewidth = op_linewidth) + 
  facet_grid(P_edge_on_cut ~ .) + 
  scale_x_continuous(breaks = pretty_breaks(n = 10)) + 
  scale_colour_manual(values = cores) + 
  theme_bw() + 
  theme(legend.position = "bottom")

g2
ggsave(plot = g2, path = wdfigs, file = 'g2.pdf', 
       width = 10, height = 10, units = "cm", device = cairo_pdf)

# TType
g3 <- ggplot(same_class) + 
  geom_density(aes(x = TType, color = misclass1), linewidth = op_linewidth) + 
  facet_grid(P_edge_on_cut ~ .) + 
  scale_x_continuous(breaks = pretty_breaks(n = 10)) + 
  scale_colour_manual(values = cores) + 
  theme_bw() + 
  theme(legend.position = "bottom")
g3
ggsave(plot = g3, path = wdfigs, file = 'g3.pdf', 
       width = 10, height = 10, units = "cm", device = cairo_pdf)

# Distância a reta
g4 <- ggplot(same_class) + 
  geom_density(aes(x = distLine_GSWLC, color = misclass1), linewidth = op_linewidth) + 
  facet_grid(P_edge_on_cut ~ .) + 
  scale_x_continuous(breaks = pretty_breaks(n = 10)) + 
  scale_colour_manual(values = cores) + 
  theme_bw() + 
  theme(legend.position = "bottom") 
g4

ggsave(plot = g4, path = wdfigs, file = 'g4.pdf', 
       width = 10, height = 10, units = "cm", device = cairo_pdf)

# vel disp
ggplot(same_class, aes(y = logvelDisp_e, x = logMstar, fill = e_disc)) + 
  geom_tile() + 
  facet_wrap(. ~ misclass1) + 
  #scale_x_continuous(breaks = pretty_breaks(n = 10)) + 
  #scale_colour_manual(values = cores) + 
  theme(legend.position = "bottom")


ggplot(faithfuld, aes(waiting, eruptions)) +
  geom_raster(aes(fill = density))











# P_S0
ggplot(same_class) + 
  geom_density(aes(x = P_S0, color = misclass1))

# P_edge_on
ggplot(same_class) + 
  geom_density(aes(x = P_edge_on, color = misclass1))

# P_cigar
ggplot(same_class) + 
  geom_density(aes(x = P_cigar, color = misclass1))

# P_S0 com TType < 0
ggplot(subset(same_class, same_class$TType < 0)) + 
  geom_density(aes(x = P_S0, color = misclass1))

# P_merg
ggplot(same_class) + 
  geom_density(aes(x = P_merg, color = misclass1))

# __B_T_r
ggplot(same_class) + 
  geom_density(aes(x = B_T_r, color = misclass1))

# dens projetada 1 vizinho
ggplot(same_class) + 
  geom_density(aes(x = logdens_proj_Neq1, color = misclass1))

# dens projetada 3 vizinho
ggplot(same_class) + 
  geom_density(aes(x = logdens_proj_Neq3, color = misclass1))

# dens projetada 5 vizinho
ggplot(same_class) + 
  geom_density(aes(x = logdens_proj_Neq5, color = misclass1))


# Preparar os dados
data(mtcars)
df <- mtcars
df$cyl <- as.factor(df$cyl)
df$am <- as.factor(df$am)

# Calcular porcentagens
library(dplyr)


# Criar o gráfico
ggplot(df_count, aes(x = bptClass, y = percentage, fill = misclass1)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal()

table(same_class$misclass1)

df_count <- same_class %>%
  group_by(flag_uv, misclass1) %>%
  summarise(count = n())

perc_corretas <- df_count$count[which(df_count$misclass1 == "Correct")]/35579
perc_FP <- df_count$count[which(df_count$misclass1 == "FP")]/4017
perc_FN <- df_count$count[which(df_count$misclass1 == "FN")]/4699

df_count$perc <- NA
df_count$perc[which(df_count$misclass1 == "Correct")] <- perc_corretas * 100
df_count$perc[which(df_count$misclass1 == "FP")] <- perc_FP * 100
df_count$perc[which(df_count$misclass1 == "FN")] <- perc_FN * 100

df_count$flag_uv <- factor(df_count$flag_uv)
levels(df_count$flag_uv) <- c("no UV", "FUV only", "NUV only", "Both")

# flag_uv
ggplot(df_count) + 
  geom_col(aes(x = flag_uv, y = perc, color = misclass1, fill = misclass1), 
           position = "dodge") + 
  facet_wrap(. ~ P_edge_on_cut) + 
  scale_x_continuous(breaks = pretty_breaks(n = 10)) + 
  scale_colour_manual(values = cores) + 
  theme(
    legend.position = "bottom")

# bptClass
ggplot(same_class) + 
  geom_bar(aes(x = bptClass, color = misclass1))