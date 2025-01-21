# SF OR Quiescent classification 

SF_Q_class <- function(logSSFR, logMstar, slope_class, intercept_class){
  
  SF <- logSSFR > (slope_class*logMstar) + intercept_class
  SF <- ifelse(SF == TRUE, "Star-forming", "Quiescent")

  return(SF)
}
