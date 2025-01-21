# dados 'clean': -------------
`%!in%` = Negate(`%in%`)
original <- read_csv("~/Work/Research/Data/SDSS_DR18_Legacy_MGS_QSO_GSWLC-X2_MGS_Simard11_Lim17_clean.csv")
# Removendo satélites que não tem no Lim17
original <- original[-which(is.na(original$igal_Lim17)),]

length(unique(original$groupID_L)) # 175894 grupos
length(which(original$centralFLAG_L)) # 171171 centrais
length(unique(original$groupID_L)) - length(which(original$centralFLAG_L)) # 4723 grupos sem centrais

centrais <- original[which(original$centralFLAG_L),]

# tabela de grupos que não tem centrais
grupos_sem_centrais <- original[which(original$groupID_L %!in% centrais$groupID_L),]
IDs_grupos_sem_centrais <- unique(grupos_sem_centrais$groupID_L) # IDs de grupos que não tem centrais

# Achar esses grupos sem centrais no Lim17
lim <- read_csv("~/Work/Research/Data/Raw data/Lim+17/SDSS(ML) galaxy - groups_M_R_c_maite.csv")
grupos_sem_centrais_LIM <- lim[which(lim$groupID_L %in% IDs_grupos_sem_centrais),]
length(unique(grupos_sem_centrais_LIM$groupID_L)) # tudo ok!

# Somente as centrais dos 'grupos sem centrais'
centrais_grupos_sem_centrais <- grupos_sem_centrais_LIM[which(grupos_sem_centrais_LIM$centralFLAG_L == TRUE),]
IDs_centrais_grupos_sem_centrais <- centrais_grupos_sem_centrais$groupID_L

write.csv(centrais_grupos_sem_centrais, 
          "/home/maitemuckler/Work/Research/Data/centrais_grupos_sem_centrais_clean.csv", row.names = F)


original_grupos_sem_centrais <- original[original$groupID_L %in% IDs_centrais_grupos_sem_centrais,]


# testando
groupID_L_1 <- original[which(original$groupID_L == 1),]
table(groupID_L_1$centralFLAG_L)




