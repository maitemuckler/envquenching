library(data.table)

input_file <- "~/Work/Research/Astronomy/Data/InputModel/inputdata_zmax0.1_Rlim2.5_Ma12.3_flag_good==1_MANGLE_logMstar_min10.5.csv"
df <- fread(input_file)

# Somente satélites:
df <- subset(df, df$type == "Satellite")
df <- df %>%
  rename(SF = SF_GSWLC)

# Deixando só colunas que importam:
variaveis <- c("SF",
               "AGN",
               "logvelDisp_e",
               "logMstar",
               "B_T_r",
               "TType",
               "conc",
               "d4000_n",
               "logMgroup",
               "logRproj_rvir",
               "logdens_n_Rproj0.5",
               "logdens_n_Rproj1",
               "logdens_proj_Neq1",
               "logdens_proj_Neq3",
               "logdens_proj_Neq5")

var_int <- c("SF",
             "AGN",
             "logvelDisp_e",
             "logMstar",
             "B_T_r",
             "TType",
             "conc",
             "d4000_n")

var_ext <- c("SF",
             "logMgroup",
             "logRproj_rvir",
             "logdens_n_Rproj0.5",
             "logdens_n_Rproj1",
             "logdens_proj_Neq1",
             "logdens_proj_Neq3",
             "logdens_proj_Neq5")

df <- df %>%
  select(all_of(variaveis))

# Somente não-agn:
df <- subset(df, df$AGN == "Non-AGN")

# Modelagem:
df$SF      <- ifelse(df$SF == "Star-forming", 1, 0)
trainIndex <- createDataPartition(df$SF, p = 0.7, list = FALSE)
train_data <- df[trainIndex, ]
test_data  <- df[-trainIndex, ]

rm(trainIndex)
rm(df)

modelo1 <- glm(SF ~ logMstar + logvelDisp_e + logRproj_rvir, family = binomial(link = "logit"), data = train_data)
