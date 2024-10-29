library(imager)
### Function: selectCoord(image)
# Allows the user to select multiple points on the provided image by clicking on it.
# Returns a list containing the coordinates of the selected points with two elements: x and y.
selectCoord <- function(image) {
  plot(image)
  coordinates <- locator(type = "p")
  return(coordinates)
}

### Function: calcParameters(solution, xmax, ymax)
# Given the solution and the image dimensions (xmax, ymax),
# this function calculates all necessary parameters for use in semla functions.
calcParameters <- function(solution, xmax, ymax) {
  # Calculate angles and restrict cosine and sine values to the range [-1, 1]
  coseno <- solution[["coseno"]]
  seno <- solution[["seno"]]
  if (coseno > 1) {coseno <- 1} else {if (coseno < -1) {coseno <- -1}}
  if (seno > 1) {seno <- 1} else {if (seno < -1) {seno <- -1}}
  
  tangente <- seno/coseno
  angulo <- atan2(seno, coseno) * (180/pi)
  
  # Calculate x and y translation parameters
  trx <- solution[["dx"]] / xmax
  try <- - solution[["dy"]] / ymax
  
  mirrorx <- solution[["mirrorx"]]
  mirrory <- solution[["mirrory"]]
  
  valores <- c(angulo, trx, try, mirrorx, mirrory)
  names(valores) <- c('angulo', 'trx', 'try', 'mirrorx', 'mirrory')
  
  # Replace any NA values with zero
  for ( i in seq_along(valores) ) {
    if (is.na(valores[[i]])) {valores[[i]] <- 0}
  }
  
  return(valores)
}

### Function: calcParametersScale(solution, xmax, ymax)
# Given the solution, xmax, and ymax parameters, this function calculates
# all necessary parameters for semla functions, including scaling.
calcParametersScale <- function(solution, xmax, ymax) {
  # Calculate angles and restrict cosine and sine values to the range [-1, 1]
  coseno <- solution[["coseno"]]
  seno <- solution[["seno"]]
  if (coseno > 1) {coseno <- 1} else {if (coseno < -1) {coseno <- -1}}
  if (seno > 1) {seno <- 1} else {if (seno < -1) {seno <- -1}}
  
  tangente <- seno/coseno
  angulo <- atan2(seno, coseno) * (180/pi)
  
  # Calculate x and y translation parameters and scaling factor
  trx <-  solution[["dx"]] / xmax
  try <- - solution[["dy"]] / ymax
  
  e <- solution[["e"]]  # Scaling factor
  
  mirrorx <- solution[["mirrorx"]]
  mirrory <- solution[["mirrory"]]
  
  valores <- c(angulo, trx, try, e, mirrorx, mirrory)
  names(valores) <- c('angulo', 'trx', 'try', 'e', 'mirrorx', 'mirrory')
  
  # Replace any NA values with zero
  for ( i in seq_along(valores) ) {
    if (is.na(valores[[i]])) {valores[[i]] <- 0}
  }
  
  return(valores)
}
