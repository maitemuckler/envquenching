# WHAN diagram

WHAN_class <- function(logN2Ha, EW_Ha_6562, EW_NII_6583){
  
  SF   <- logN2Ha < -0.4 & EW_Ha_6562 > 3
  sAGN <- logN2Ha > -0.4 & EW_Ha_6562 > 6
  wAGN <- logN2Ha > -0.4 & EW_Ha_6562 > 3 & EW_Ha_6562 < 6
  RG   <- EW_Ha_6562 < 3
  PG   <- EW_Ha_6562 < 0.5 & EW_NII_6583 < 0.5
  
  whan <- NA
  whan <- ifelse(SF == TRUE, 'SF', whan)
  whan <- ifelse(sAGN == TRUE, 'sAGN', whan)
  whan <- ifelse(wAGN == TRUE, 'wAGN', whan)
  whan <- ifelse(RG == TRUE, 'RG', whan)
  whan <- ifelse(PG == TRUE, 'PG', whan)
  
  return(whan)
}