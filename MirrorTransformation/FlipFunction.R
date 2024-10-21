# Funcion inversion en el eje Y (newrow = 77 - row)
inversionEjeY <- function(muestra1) {
  # se guarda el dataframe de las coordenadas de la muestra 1
  df1 <- muestra1@images$slice1@coordinates
  # ordenar el data.frame
  df1ord <- df1[order(df1$row, df1$col),]
  
  # se gurdan los imagenrow/imagecol de cada row/col
  rowRelation <- list()
  colRelation <- list()
  
  for ( barcode in rownames(df1ord) ) { 
    rowB <- as.character(df1ord[barcode,"row"])
    imagerowB <- df1ord[barcode,"imagerow"]
    
    colB <- as.character(df1ord[barcode,"col"])
    imagecolB <- df1ord[barcode,"imagecol"]
    
    if ( !(imagerowB %in% rowRelation) ) {
      rowRelation[[rowB]] <- imagerowB
    } #end_if
    
    if ( !(imagecolB %in% colRelation) ) {
      colRelation[[colB]] <- imagecolB
    } #end_if
  } #end_for
  
  # se crea df2 = df1
  df2 <- df1
  for ( barcode in rownames(df1) ) { 
    # se modifica sus valores de row y de imagerow 
    rowB <- df1[barcode,"row"]
    colB <- df1[barcode,"col"]
    newrowB <- 77 - rowB
    
    df2[barcode,"row"] <- newrowB
    df2[barcode,"imagerow"] <- rowRelation[[as.character(newrowB)]]
  } #end_for
  
  # se crea la nueva muestra modificada
  muestra2 <- muestra1
  nameNew <- c()
  for ( i in 1:length(muestra1@meta.data$name) ) { nameNew <- c(nameNew, paste0(muestra1@meta.data$name[i], "_Xinvert")) }
  muestra2@meta.data$name <- nameNew
  muestra2@images$slice1@coordinates <- df2
  
  imagen_volteada <- array(NA, dim = dim(muestra1@images$slice1@image))
  for (i in 1:dim(muestra1@images$slice1@image)[3]) {
    imagen_volteada[,,i] <- muestra1@images$slice1@image[nrow(muestra1@images$slice1@image):1, , i]
  }
  
  muestra2@images$slice1@image <- imagen_volteada
  
  return(muestra2)
} #end_function


# Funcion inversion en el eje X (newcol = 127 - col)
inversionEjeX <- function(muestra1) {
  # se guarda el dataframe de las coordenadas de la muestra 1
  df1 <- muestra1@images$slice1@coordinates
  # ordenar el data.frame
  df1ord <- df1[order(df1$row, df1$col),]
  
  # se gurdan los imagenrow/imagecol de cada row/col
  rowRelation <- list()
  colRelation <- list()
  
  for ( barcode in rownames(df1ord) ) { 
    rowB <- as.character(df1ord[barcode,"row"])
    imagerowB <- df1ord[barcode,"imagerow"]
    
    colB <- as.character(df1ord[barcode,"col"])
    imagecolB <- df1ord[barcode,"imagecol"]
    
    if ( !(imagerowB %in% rowRelation) ) {
      rowRelation[[rowB]] <- imagerowB
    } #end_if
    
    if ( !(imagecolB %in% colRelation) ) {
      colRelation[[colB]] <- imagecolB
    } #end_if
  } #end_for
  
  # se crea df2 = df1
  df2 <- df1
  for ( barcode in rownames(df1) ) { 
    # se modifica sus valores de col y de imagecol 
    rowB <- df1[barcode,"row"]
    colB <- df1[barcode,"col"]
    newcolB <- 127 - colB
    
    df2[barcode,"col"] <- newcolB
    df2[barcode,"imagecol"] <- colRelation[[as.character(newcolB)]]
  } #end_for
  
  # se crea la nueva muestra modificada
  muestra2 <- muestra1
  nameNew <- c()
  for ( i in 1:length(muestra1@meta.data$name) ) { nameNew <- c(nameNew, paste0(muestra1@meta.data$name[i], "_Yinvert")) }
  muestra2@meta.data$name <- nameNew
  muestra2@images$slice1@coordinates <- df2
  
  imagen_volteada <- array(NA, dim = dim(muestra1@images$slice1@image))
  for (i in 1:dim(muestra1@images$slice1@image)[3]) {
    imagen_volteada[,,i] <- muestra1@images$slice1@image[, ncol(muestra1@images$slice1@image):1, i]
  }
  
  muestra2@images$slice1@image <- imagen_volteada
  
  return(muestra2)
} #end_function


# Funcion inversion en el eje X (newrow = 77 - row) e Y (newcol = 127 - col)
inversionEjeXY <- function(muestra1) {
  # se guarda el dataframe de las coordenadas de la muestra 1
  df1 <- muestra1@images$slice1@coordinates
  # ordenar el data.frame
  df1ord <- df1[order(df1$row, df1$col),]
  
  # se gurdan los imagenrow/imagecol de cada row/col
  rowRelation <- list()
  colRelation <- list()
  
  for ( barcode in rownames(df1ord) ) { 
    rowB <- as.character(df1ord[barcode,"row"])
    imagerowB <- df1ord[barcode,"imagerow"]
    
    colB <- as.character(df1ord[barcode,"col"])
    imagecolB <- df1ord[barcode,"imagecol"]
    
    if ( !(imagerowB %in% rowRelation) ) {
      rowRelation[[rowB]] <- imagerowB
    } #end_if
    
    if ( !(imagecolB %in% colRelation) ) {
      colRelation[[colB]] <- imagecolB
    } #end_if
  } #end_for
  
  # se crea df2 = df1
  df2 <- df1
  for ( barcode in rownames(df1) ) { 
    # se modifica sus valores de row y col y de image row y imagecol 
    rowB <- df1[barcode,"row"]
    colB <- df1[barcode,"col"]
    newrowB <- 77 - rowB
    newcolB <- 127 - colB
    
    df2[barcode,"row"] <- newrowB
    df2[barcode,"imagerow"] <- rowRelation[[as.character(newrowB)]]
    
    df2[barcode,"col"] <- newcolB
    df2[barcode,"imagecol"] <- colRelation[[as.character(newcolB)]]
  } #end_for
  
  # se crea la nueva muestra modificada
  muestra2 <- muestra1
  nameNew <- c()
  for ( i in 1:length(muestra1@meta.data$name) ) { nameNew <- c(nameNew, paste0(muestra1@meta.data$name[i], "_XYinvert")) }
  muestra2@meta.data$name <- nameNew
  muestra2@images$slice1@coordinates <- df2
  
  imagen_volteada <- array(NA, dim = dim(muestra1@images$slice1@image))
  for (i in 1:dim(muestra1@images$slice1@image)[3]) {
    imagen_volteada[,,i] <- muestra1@images$slice1@image[nrow(muestra1@images$slice1@image):1, ncol(muestra1@images$slice1@image):1, i]
  }
  
  muestra2@images$slice1@image <- imagen_volteada
  
  return(muestra2)
} #end_function


