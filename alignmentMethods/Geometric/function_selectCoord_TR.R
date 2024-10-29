### Function: selectCoord(image)
# Allows the user to select multiple points on the provided image by clicking on it.
# Returns a list containing the coordinates of the selected points with two elements: x and y.
selectCoord <- function(image) {
  library(imager)
  plot(image)
  coordinates <- locator(type = "p")
  return(coordinates)
}

### Function: solveCoord(coord1, coord2, xmax, ymax)
# Given two sets of coordinates (coord1 and coord2), and the maximum x and y values (xmax, ymax),
# this function calculates cos(angle), sin(angle), dx, and dy by solving a system of equations.
solveCoord <- function(coord1, coord2, xmax, ymax) {
  cx <- xmax / 2
  cy <- ymax / 2
  
  # Adjust coordinates to be centered based on xmax and ymax
  coord2MOD <- coord2
  for (j in seq_along(coord2$x)) {
    coord2MOD$x[[j]] <- coord2$x[[j]] - cx
    coord2MOD$y[[j]] <- coord2$y[[j]] - cy
  }
  
  coord1MOD <- coord1
  for (j in seq_along(coord2$x)) {
    coord1MOD$x[[j]] <- coord1$x[[j]] - cx
    coord1MOD$y[[j]] <- coord1$y[[j]] - cy
  }
  
  # Optimization of an overdetermined equation system
  fn <- function(x) {
    y <- numeric(0)
    for (i in seq_along(coord1$x)) {
      y[length(y) + 1] <- coord2MOD$x[[i]]*x[1] + coord2MOD$y[[i]]*x[2] + 1*x[3] + 0*x[4] - coord1MOD$x[[i]]
      y[length(y) + 1] <- coord2MOD$y[[i]]*x[1] - coord2MOD$x[[i]]*x[2] + 0*x[3] + 1*x[4] - coord1MOD$y[[i]]
    }
    y[length(y) + 1] <- (x[1])^2 + (x[2])^2 - 1  
    return(sum(y^2))
  }
  
  xstart <- c(cos(pi/4), sin(pi/4), 2, -2)
  solution <- optim(xstart, fn)
  solution <- solution$par
  names(solution) <- c('coseno', 'seno', 'dx', 'dy')
  
  return(solution)
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
  
  tangente <- seno / coseno
  angulo <- atan2(seno, coseno) * (180/pi)
  
  # Calculate x and y translation parameters
  trx <- solution[["dx"]] / xmax
  try <- - solution[["dy"]] / ymax
  
  mirrorx <- solution[["mirrorx"]]
  mirrory <- solution[["mirrory"]]
  
  valores <- c(angulo, trx, try, mirrorx, mirrory)
  names(valores) <- c('angulo', 'trx', 'try', 'mirrorx', 'mirrory')
  
  # Replace any NA values with zero
  for (i in seq_along(valores)) {
    if (is.na(valores[[i]])) {valores[[i]] <- 0}
  }
  
  return(valores)
}

### Function: calcNewCoord(coord, sol, xmax, ymax)
# Given the original coordinates (coord), calculated transformation parameters (sol),
# and image dimensions (xmax, ymax), this function calculates the new transformed coordinates.
calcNewCoord <- function(coord, sol, xmax, ymax) {
  cx <- xmax / 2
  cy <- ymax / 2
  
  coordMOD <- coord
  
  # Apply mirror transformations based on solution parameters
  if (sol[["mirrorx"]] == 1) {
    for (j in seq_along(coordMOD$x)) {
      coordMOD$x[[j]] <- xmax - coordMOD$x[[j]]
    }
  }
  
  if (sol[["mirrory"]] == 1) {
    for (j in seq_along(coordMOD$y)) {
      coordMOD$y[[j]] <- ymax - coordMOD$y[[j]]
    }
  }
  
  # Center coordinates
  for (j in seq_along(coordMOD$x)) {
    coordMOD$x[[j]] <- coordMOD$x[[j]] - cx
    coordMOD$y[[j]] <- coordMOD$y[[j]] - cy
  }
  
  # Apply rotation and translation transformations
  coordNEW <- coordMOD
  for (j in seq_along(coordMOD$x)) {
    coordNEW$x[[j]] <- round(coordMOD$x[[j]]*sol[["coseno"]] + coordMOD$y[[j]]*sol[["seno"]] + sol[["dx"]] + cx)
    coordNEW$y[[j]] <- round(coordMOD$y[[j]]*sol[["coseno"]] - coordMOD$x[[j]]*sol[["seno"]] + sol[["dy"]] + cy)
  }
  
  return(coordNEW)
}

