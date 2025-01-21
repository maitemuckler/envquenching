# Unindo tabelas de vizinhos e de densidade

## Carregar pacotes ----
library(readr)
library(dplyr)
library(ggplot2)
library(scales)
library(data.table)

### Definir diretórios ----
wdmain       <- "~/Work/Research/Astronomy/"
wdproject    <- paste0(wdmain, "Projects/environmental-quenching/")
wdcode       <- paste0(wdproject, "Scripts/")
wddata       <- paste0(wdmain, "Data/")

## Carregar input data ----
input_file <- "SDSS_DR18_Legacy_MGS_QSO"

# Abertura fixa ----
ap_05 <- read_csv(paste0(wddata, paste0('LocalDensity/DensityVolume/', input_file, '_density_Rproj0.5_Dv3000kms_Mabs-20.4.csv')))
ap_10 <- read_csv(paste0(wddata, paste0('LocalDensity/DensityVolume/', input_file, '_density_Rproj1.0_Dv3000kms_Mabs-20.4.csv')))

ap <- merge(ap_05, ap_10, by = "igal")
rm(ap_05, ap_10)

# Vizinhos ----
neighbors1 <- read_csv(paste0(wddata, paste0('LocalDensity/NthSpecGalaxy/', input_file, '_Nth_spec_galaxy_Neq1_Dv3000kms_Mabs-20.4.csv')))
neighbors3 <- read_csv(paste0(wddata, paste0('LocalDensity/NthSpecGalaxy/', input_file, '_Nth_spec_galaxy_Neq3_Dv3000kms_Mabs-20.4.csv')))
neighbors5 <- read_csv(paste0(wddata, paste0('LocalDensity/NthSpecGalaxy/', input_file, '_Nth_spec_galaxy_Neq5_Dv3000kms_Mabs-20.4.csv')))

neighbors <- merge(neighbors1, neighbors3, by = "igal")
neighbors <- merge(neighbors, neighbors5, by = "igal")
rm(neighbors1, neighbors3, neighbors5)

neighbors$logdens_proj_Neq1 <- log10(neighbors$dens_proj_Neq1)
neighbors$logdens_proj_Neq3 <- log10(neighbors$dens_proj_Neq3)
neighbors$logdens_proj_Neq5 <- log10(neighbors$dens_proj_Neq5)

# Unindo abertura fixa e vizinhos ----

local_density <- merge(ap, neighbors, by = "igal")
rm(ap, neighbors)

# Lendo tabela input e unindo com densidade local ----

inputdata              <- fread(paste0(wddata, "Raw data/SDSS/", input_file, ".csv"))
inputdata_localDensity <- merge(inputdata, local_density, by = 'igal') # (108 colunas)
rm(inputdata, local_density)

output_file <- paste0(input_file, "_localDensity.csv")
write.csv(inputdata_localDensity, paste0(wddata, "Raw data/SDSS/", output_file), row.names = F)

