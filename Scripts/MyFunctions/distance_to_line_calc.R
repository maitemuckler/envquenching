distance_to_line <- function(x, y, slope, intercept) {
  A <- slope
  B <- -1
  C <- intercept
  numerator <- A * x + B * y + C
  denominator <- sqrt(A^2 + B^2)
  distance <- numerator / denominator
  return(distance)
}