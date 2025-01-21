# List the functions in this file
if(!exists('list_functions_interp')){
  print('Defining functions from interp.R:')
  print('-------------------------------------')
  print(' 1. interp(xin, xout, yin)')
  print(' 2. interp2(xin, xout, yin)')

  list_functions_interp = F
}

##################################################################################################
# FUNCTION interp
##################################################################################################
# PURPOSE:
#    
# CALLING SEQUENCE:
#    
# INPUTS:
#    
# PARAMETERS:
#    
# OUTPUT:
#    
# REQUIRED SCRIPTS:
#   
##################################################################################################  
interp <- function(xin, xout, yin){
  Ninterpol = length(xout)
  yout <- matrix(nrow = Ninterpol)
  
  for(k in 1:Ninterpol){
    t = xin[xin < xout[k]]
    tind = length(t)
    
    if(tind <= 1){tind = 2}
    if(tind >= length(xin)){tind = length(xin) - 1}
    t1 = xin[tind - 1]
    t2 = xin[tind]
    t3 = xin[tind + 1]
    tx = xout[k]
    
    A = (tx - t1) / (t3 - t1)
    B = (tx - t2) / (t3 - t2)
    C = (tx - t3) / (t2 - t1)
    D = (tx - t1) / (t3 - t2)
    E = (tx - t2) / (t3 - t1)
    
    G = (tx - t2) / (t3 - t2)
    H = (tx - t1) / (t2 - t1)
    
    ###########################################
    if(t1 != t2 & t2 != t3){
      yout[k] = (yin[tind+1] * A * B - yin[tind] * D * C + yin[tind-1] * E * C)
    }
    if(t1 == t2){
      yout[k] = (yin[tind+1] - yin[tind]) * G + yin[tind] 
    }
    if(t2 == t3){
      yout[k] = (yin[tind] - yin[tind-1]) * H + yin[tind-1] 
    }
  } ####### FIM DA INTERPOLAÇÃO
  
  return(yout)
}

##################################################################################################
# FUNCTION interp2
##################################################################################################
# PURPOSE:
#    
# CALLING SEQUENCE:
#    
# INPUTS:
#    
# PARAMETERS:
#    
# OUTPUT:
#    
# REQUIRED SCRIPTS:
#   
##################################################################################################  
interp2 <- function(xin, xout, yin){
  Ninterpol = length(xout)
  yout <- matrix(nrow = Ninterpol)
  
  for(k in 1:Ninterpol){
    t = xin[xin < xout[k]]
    tind = length(t)
    
    if(tind < 1){tind = 1}
    if(tind >= length(xin)){tind = length(xin) - 1}
    t2 = xin[tind]
    t3 = xin[tind + 1]
    tx = xout[k]

    G = (tx - t2) / (t3 - t2)
  
    ###########################################
    yout[k] = (yin[tind+1] - yin[tind]) * G + yin[tind] 
    
  } ####### FIM DA INTERPOLAÇÃO
  
  return(yout)
}

