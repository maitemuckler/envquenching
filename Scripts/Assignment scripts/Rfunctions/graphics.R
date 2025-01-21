# List the functions in this file
if(!exists('list_functions_graphics')){
  print('Defining functions from graphics.R:')
  print('-------------------------------------')
  print(' 1. polygon_vectors(x, ylo, yup)')
  print(' 2. superpose.eb(x, y, ebl, ebu, ebl, length, col, lwd)')

  list_functions_graphics = F
}

##################################################################################################
# FUNCTION poly_vectors
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
polygon_vectors <- function(x, ylo, yup){
  NN = length(x)
  x_pol = vector(length = 2 * NN)
  y_pol = vector(length = 2 * NN)

  for(i in 1:NN){
    x_pol[i] = x[i]
    x_pol[(i+NN)] = x[(NN-(i-1))]

    y_pol[i] = yup[i] 
    y_pol[(i+NN)] = ylo[(NN-(i-1))] 
  }
  
  return(data.frame(x = x_pol, y = y_pol))
}

##################################################################################################
# FUNCTION superpose.eb
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
superpose.eb <- function(x, y, ebl, ebu, length, col, lwd){
  if(missing(ebu)){ebu = ebl}
  if(missing(col)){col = 'black'}
  if(missing(lwd)){lwd = 1}
  if(missing(length)){
    if(length(x) > 1){
      length = (max(x) - min(x)) / 100
    }else{
      length = 0.02
    }
  }
  arrows(x, y + ebl, x, y - ebu, angle = 90, code = 3, length = length, col = col, lwd = lwd)
}

