library(recipes)
library(RColorBrewer)
library(wesanderson)
library(scico)
library(ggpubr)
library(ggthemes)


# Contagem de grupos, satélites, etc ----

# Satélites vs central
table(df$type)
table(df$type)/nrow(df)*100

# Grupos
length(unique(df$groupID))

# Centrais
length(which(df$Rproj_rvir == 0))


