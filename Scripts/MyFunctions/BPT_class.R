# BPT diagram

# seleciona quem não é AGN
BPT_class <- function(logO3Hb, logN2Ha, aa, bb){
  
  AGN <- logO3Hb <= 0.61 / (logN2Ha - aa) + bb &
         logO3Hb > -3 & logO3Hb < 2 &
         logN2Ha > -3 & logN2Ha < 0.4
  
  AGN[which(is.na(AGN))] <- TRUE
  
  AGN[which(AGN == TRUE)]  <- "Non-AGN"
  AGN[which(AGN == FALSE)] <- "AGN"

  return(AGN)
  
}
