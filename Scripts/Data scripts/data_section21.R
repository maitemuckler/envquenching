# log10(Mstar/Mhalo): 
df$Mstar_Mhalo <- df$logMstar - df$logMgroup

## Criar colunas para variáveis a respeito das centrais ----



# Colunas que quero salvar a informação da central 
colsCentral <- c("groupID", 
                 "logMstar", 
                 "logSFR_SED", 
                 "logvelDisp_e", 
                 "vlos_vvir", 
                 "conc", 
                 "d4000_n",
                 "logSigma_SFR", 
                 "distLine_GSWLC", 
                 "SF_GSWLC", 
                 "AGN",
                 "P_disk", 
                 "P_edge_on", 
                 "P_bar_GZ2", 
                 "P_bar_Nair10", 
                 "P_merg", 
                 "P_bulge", 
                 "B_T_r", 
                 "e", 
                 "P_cigar", 
                 "TType", 
                 "P_S0", 
                 "TType_label")

colsCentral_names <- colsCentral
colsCentral_names[-1] <- paste0(colsCentral[-1], "_central")

aux <- df %>%
  filter(type == "Central") %>%
  select(all_of(colsCentral))

colnames(aux) <- colsCentral_names

df <- merge(df, aux, by = "groupID", all.x = T)

rowIndex <- which(colnames(df) %in% colsCentral_names)
rowIndex <- rowIndex[-1]
nrow(df[rowSums(is.na(df[,..rowIndex])) == ncol(df[,..rowIndex]), ..rowIndex]) # Sem informação de central

rm(aux, df_assign, df_clean)

# Variáveis que quero transformar para factor:
colunas_para_converter <- c("Class",
                            "bptClass",
                            "flag_good",
                            "flag_sed",
                            "flag_uv",
                            "flag_midir",
                            "type",
                            "SF_GSWLC",
                            "SF_MPAJHU",
                            "TType_label",
                            "AGN",
                            "groupID", 
                            "SF_GSWLC_central", 
                            "AGN_central",
                            "TType_label_central")


## Verificando valores NA e inf ----

names(df)[sapply(df, function(x) any(is.na(x)))]
names(df)[sapply(df, function(x) any(is.infinite(x)))]

# logRproj_rvir        <- infinito porque era central, então Rproj = 0 -> logRproj = Inf.
# logO3Hb e logN2Ha    <- alguma divisão por zero.
# logvelDisp_e_central <- a central deve ter velDisp_e = 0 -> logvelDisp_e = Inf.
# Para logvelDisp_e_central vou susbtituir o Inf por NA.

# Central
df$logvelDisp_e_central[which(is.infinite(df$logvelDisp_e_central))] <- NA