# SF, Quiescent OR GV classification 

SF_Q_GV_class <- function(logSSFR, logMstar, slope_class, intercept_class, 
                          sigma_GV_up, sigma_GV_down){
  
  # Equation 
  m <- (slope_class*logMstar) + intercept_class
  
  star_forming <- logSSFR > m + (sigma_GV_up*sd(m))
  quiescent    <- logSSFR < m - (sigma_GV_down*sd(m))
  
  SF_Q_GV <- NA
  SF_Q_GV <- ifelse(star_forming == TRUE, 'Star-forming', SF_Q_GV)
  SF_Q_GV <- ifelse(quiescent == TRUE, 'Quiescent', SF_Q_GV)
  SF_Q_GV[which(is.na(SF_Q_GV))] <- 'Green Valley'
  
  return(SF_Q_GV)
}

