### function selectCoord(image)
# Allows the user to select multiple points on the provided image by clicking on it.
# Returns a list containing the coordinates of the selected points with two elements: x and y.
selectCoord <- function(image) {
  library(imager)
  plot(image)
  coordinates <- locator(type = "p")
  return(coordinates)
}

# Load required libraries
library(reticulate)
library(semla)
library(tibble)

# Define directories for Seurat and semla objects
carpetaData <- "/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/02.Seurat/Yadav2023_adult_human_spinal_cord/ImagenesYadavMergeSlices/"
saveDir <- "/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/05.ImageJ"

# Load a custom function for alignment evaluation
source("/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/04.ImageAlingment/EvalAlign/function_evalAlign.R")

# Load the merged Seurat object 
paciente19.merge <- readRDS(file = paste0(carpetaData, "/Paciente19_merge.rds"))

# Update Seurat object for compatibility with the semla package
pacientes.semla <- UpdateSeuratForSemla(paciente19.merge)

# Load images into the Seurat object with specified image height
pacientes.semla <- LoadImages(pacientes.semla, image_height = 600)

# Plot loaded images for visual inspection
ImagePlot(pacientes.semla)

# Set up lists for alignment process and coordinate retrieval
aligned <- list(pacientes.semla@tools$Staffli@rasterlists$raw[[1]])
original <- list(pacientes.semla@tools$Staffli@rasterlists$raw[[1]])
listaSoluciones <- list(NA)
listaTransforms <- list(NA)

# Extract column and row coordinates for sample ID 1 from Seurat object
CoordinatesSemlaCol <- list(GetCoordinates(pacientes.semla)[which(GetCoordinates(pacientes.semla)["sampleID"] == 1), 2])
CoordinatesSemlaRow <- list(GetCoordinates(pacientes.semla)[which(GetCoordinates(pacientes.semla)["sampleID"] == 1), 3])


"
This is the point where you would save the images from the RDS file into a source_dir folder, 
create a target_dir folder, and then use ImageJ with the Register Virtual Stack Slices plugin for alignment.
"

# Define source and target directories for ImageJ alignment
source_dir <- "/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/05.ImageJ/ImagenesAlinearFiji/paciente19"
target_dir <- "/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/05.ImageJ/ImagenesAlinearFiji/paciente19aligned"
#ref_name <- "original"
ref_name <- "tissue_lowres_image_1"

# Instructions for folder setup before alignment
cat(paste("You should have a folder called: ", source_dir,
          "\nAnd a folder called: ", target_dir,
          "\nIn the first one, you should have the images from de rds document."))

# Provide instructions for the manual alignment process in ImageJ
cat(paste("It is time to open ImageJ and do the alignment of the images.",
          "\nFirst, go to File>Open and select the different images to align.",
          "The images should have only one color channel.",
          "To do this you should go to Image>Color>Split Channels and then Image>Color>Merge Channels unselecting create composite.",
          "\nTo do the align go to Plugins>Registration>Register Virtual Stack Slices.", 
          "Remember you should save the transforms.",
          "\nThe transformation parameters or Transform files should be stored in the same folder as the result images.",
          "\nAnd finally you have to select the reference image."))

# Wait for user confirmation before proceeding
answer <- "no"
while (answer != 'yes') {
  cat(paste0("When you have all this completed, write: 'yes': "))
  answer <- readline()
}
if (answer == 'yes') {
cat(paste0("The alignment is done and the transformation parameters are saved."))
}

# Extract translation and rotation parameters from XML files in the target directory
filesTarget <- list.files(target_dir)
files_xml <- filesTarget[grep("\\.xml$", filesTarget)]

Nimage <- 1
for (xml in files_xml) {
  # Skip the reference image file and process each alignment XML file
  if (!startsWith(xml, ref_name)) {
    Nimage <- Nimage + 1
    lines <- readLines(paste0(target_dir, "/", xml))
    lines <- lines[startsWith(lines, "\t<iict_transform")]
    lines <- sapply(lines, function(line) substr(line, start = nchar("\\t<iict_transform"), stop = nchar(line) - 3))
    lines <- sapply(lines, function(line) strsplit(line, "\""))
    
    # Extract rotation and translation parameters based on transform class type
    for (line in lines) {
      classN <- which(startsWith(line, " class="))
      if (endsWith(line[[classN + 1]], "transform.RigidModel2D")) { 
        dataN <- which(startsWith(line, " data="))
        transformsR <- line[[dataN + 1]]
        transformsR <- strsplit(transformsR, " ")
      }
      if (endsWith(line[[classN + 1]], "transform.TranslationModel2D")) { 
        dataN <- which(startsWith(line, " data="))
        transformsT <- line[[dataN + 1]]
        transformsT <- strsplit(transformsT, " ")
      }
    }
    # Save transformation parameters: rotation angle, dx, and dy
    transformsParams <- list()
    transformsParams["angle"] <- as.numeric(transformsR[[1]][1])  # radians
    transformsParams["dx"] <- as.numeric(transformsR[[1]][2]) + as.numeric(transformsT[[1]][1])  # without normalization (-1,1)
    transformsParams["dy"] <- as.numeric(transformsR[[1]][3]) + as.numeric(transformsT[[1]][2])  # without normalization (-1,1)
    listaSoluciones[[Nimage]] <- transformsParams
  }
}

# Apply transformations for each image in the Seurat object
for (i in 2:length(pacientes.semla@tools$Staffli@rasterlists$raw)) {
  # Set image dimensions based on minimum height and width across images
  xmax <- min(nrow(pacientes.semla@tools$Staffli@rasterlists$raw[[1]]),
              nrow(pacientes.semla@tools$Staffli@rasterlists$raw[[i]]))
  ymax <- min(ncol(pacientes.semla@tools$Staffli@rasterlists$raw[[1]]),
              ncol(pacientes.semla@tools$Staffli@rasterlists$raw[[i]]))
  
  # Retrieve and convert transformation parameters
  solucion <- listaSoluciones[[i]]
  valores <- solucion
  valores$angle <- solucion$angle * (180/pi)  # Convert angle to degrees
  
  # Adjust dx and dy based on image center for transformation normalization
  valores$dx <-  (solucion$dx - (xmax/2 - xmax/2 * cos(solucion$angle) + ymax/2 * sin(solucion$angle))) / xmax
  valores$dy <-  (solucion$dy - (ymax/2 - ymax/2 * cos(solucion$angle) - xmax/2 * sin(solucion$angle))) / ymax
  
  alltransforms <- list()
  original[[i]] <- pacientes.semla@tools$Staffli@rasterlists$raw[[i]]
  
  
  # Apply rotation if the angle is non-zero
  if (valores$angle != 0 && valores$angle != 360) {
    transforms_angle <- generate_rigid_transform(sampleID = i, 
                                                 angle = valores$angle)
    pacientes.semla <- RigidTransformImages(pacientes.semla, transforms = transforms_angle)
    alltransforms[[2]] <- transforms_angle
    pacientes.semla@tools$Staffli@rasterlists$raw[[i]] <- pacientes.semla@tools$Staffli@rasterlists$transformed[[i]]
    
    # Update meta_data with new transformations
    pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 2] <- pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 4]
    pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 3] <- pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 5]
  }
  
  # Apply translation if dx or dy is non-zero
  if (valores$dx != 0 || valores$dy != 0) {
    transforms_trans <- generate_rigid_transform(sampleID = i, 
                                                 tr_x = valores$dx, 
                                                 tr_y = valores$dy)
    pacientes.semla <- RigidTransformImages(pacientes.semla, transforms = transforms_trans)
    alltransforms[[3]] <- transforms_trans
    pacientes.semla@tools$Staffli@rasterlists$raw[[i]] <- pacientes.semla@tools$Staffli@rasterlists$transformed[[i]]
    
    # Update meta_data with new transformations
    pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 2] <- pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 4]
    pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 3] <- pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 5]
  }
  
  # Store aligned images and transformation lists
  aligned[[i]] <- pacientes.semla@tools$Staffli@rasterlists$transformed[[i]]
  pacientes.semla@tools$Staffli@rasterlists$raw[[i]] <- original[[i]]
  listaTransforms[[i]] <- alltransforms
  
  # Save transformed coordinates for each sample
  CoordinatesSemlaCol[[i]] <- GetCoordinates(pacientes.semla)[which(GetCoordinates(pacientes.semla)[6] == i), 4]
  CoordinatesSemlaRow[[i]] <- GetCoordinates(pacientes.semla)[which(GetCoordinates(pacientes.semla)[6] == i), 5]
}

# Update transformed image list in Seurat object
pacientes.semla@tools$Staffli@rasterlists[["transformed"]] <- aligned

# Aggregate transformed column and row coordinates
CoordCol <- integer()
CoordRow <- integer()
for (i in seq_along(CoordinatesSemlaCol)) {
  CoordCol <- c(CoordCol, CoordinatesSemlaCol[[i]][[1]])
  CoordRow <- c(CoordRow, CoordinatesSemlaRow[[i]][[1]])
}
pacientes.semla@tools$Staffli@meta_data$pxl_col_in_fullres_transformed <- CoordCol
pacientes.semla@tools$Staffli@meta_data$pxl_row_in_fullres_transformed <- CoordRow

# Plot the original and transformed images
ImagePlot(pacientes.semla)
title("Original H&E image")
png(paste0(saveDir, "/results/patient19_merge_imageJ.png"))
ImagePlot(pacientes.semla, image_use = "transformed")
title("Transformed H&E image", line = -12.8)
dev.off()


# Evaluate the alignment
# Initialize a list to store coordinates
listaCoordenadas <- list()

# Loop through the raw raster lists of the 'pacientes.semla' object to extract coordinates
for (i in seq_along(pacientes.semla@tools$Staffli@rasterlists$raw)) {
  coordenadas <- selectCoord(pacientes.semla@tools$Staffli@rasterlists$raw[[i]])
  
  # Round the coordinates to the nearest integer
  for (j in seq_along(coordenadas)) {
    coordenadas[[j]] <- round(coordenadas[[j]])
  }
  
  # Store the rounded coordinates in the list
  listaCoordenadas[[i]] <- coordenadas
}

# Create a new list for transformed images' coordinates
listaCoordenadasNEW <- listaCoordenadas

# Extract coordinates for transformed images, similar to the raw images
for (i in 2:length(pacientes.semla@tools$Staffli@rasterlists$transformed)) {
  coordenadas <- selectCoord(pacientes.semla@tools$Staffli@rasterlists$transformed[[i]])
  
  # Round the coordinates
  for (j in seq_along(coordenadas)) {
    coordenadas[[j]] <- round(coordenadas[[j]])
  }
  
  # Store the rounded coordinates in the new list
  listaCoordenadasNEW[[i]] <- coordenadas
}

# Initialize an evaluation list
Evaluation <- list()
# Define coordinates for comparison
x1 <- listaCoordenadasNEW[[1]]$x[[5]]
y1 <- listaCoordenadasNEW[[1]]$y[[5]]

# Create a list to store control alignment parameters
control <- list()

# Loop through the new coordinates and evaluate the alignment
for (i in 1:length(listaCoordenadasNEW)) {
  control[[i]] <- controlAlign(pacientes.semla@tools[["Staffli"]]@rasterlists[["transformed"]][[i]])
  print(control)
}

# Store the control evaluations in the Evaluation list
Evaluation$control <- control

# Evaluate raw images (same region, same point)
for ( i in 1:length(listaCoordenadas)) {
  if (i != length(listaCoordenadas)) {
    for ( j in i:length(listaCoordenadas))  {
      if (i != j) {
        x1 <- listaCoordenadas[[i]]$x[[5]]
        y1 <- listaCoordenadas[[i]]$y[[5]]
        print(paste0("Evaluation of the alignment between images ", as.character(i), " and ", as.character(j), " comparing the same region without selecting a different point."))
        
        # Evaluate the alignment and store parameters
        parameters <- 
          evalAlign(
            pacientes.semla@tools[["Staffli"]]@rasterlists[["raw"]][[i]][(x1-200):(x1+200),(y1-200):(y1+200)],
            pacientes.semla@tools[["Staffli"]]@rasterlists[["raw"]][[j]][(x1-200):(x1+200),(y1-200):(y1+200)], 
            listaCoordenadas, c(i,j))
        print(parameters)
        imagescompare <- paste0("comparing_", as.character(i), '_', as.character(j))
        Evaluation$original$sameRegion_samePoint[[imagescompare]] <- parameters
        
      }}}}

# Evaluate raw images (common region, different point)
for ( i in 1:length(listaCoordenadas)) {
  if (i != length(listaCoordenadas)) {
    for ( j in i:length(listaCoordenadas))  {
      if (i != j) {
        x1 <- listaCoordenadas[[i]]$x[[5]]
        y1 <- listaCoordenadas[[i]]$y[[5]]
        x2 <- listaCoordenadas[[j]]$x[[5]]
        y2 <- listaCoordenadas[[j]]$y[[5]]
        print(paste0("Evaluation of the alignment between images ", as.character(i), " and ", as.character(j), " selecting an area from a reference point in each image."))
        
        # Evaluate the alignment and store parameters
        parameters <-
          evalAlign(
            pacientes.semla@tools[["Staffli"]]@rasterlists[["raw"]][[i]][(x1-200):(x1+200),(y1-200):(y1+200)],
            pacientes.semla@tools[["Staffli"]]@rasterlists[["raw"]][[j]][(x2-200):(x2+200),(y2-200):(y2+200)],
            listaCoordenadas, c(i,j))
        print(parameters)
        imagescompare <- paste0("comparing_", as.character(i), '_', as.character(j))
        Evaluation$original$commonRegion_differentPoint[[imagescompare]] <- parameters
      }}}}

# Evaluate transformed images (same region, same point)
for ( i in 1:length(listaCoordenadasNEW)) {
  if (i != length(listaCoordenadasNEW)) {
    for ( j in i:length(listaCoordenadasNEW))  {
      if (i != j) {
        x1 <- listaCoordenadasNEW[[i]]$x[[5]]
        y1 <- listaCoordenadasNEW[[i]]$y[[5]]
        print(paste0("Evaluation of the alignment between images ", as.character(i), " and ", as.character(j), " comparing the same region without selecting a different point."))
        
        # Evaluate the alignment and store parameters
        parameters <- 
          evalAlign(
            pacientes.semla@tools[["Staffli"]]@rasterlists[["transformed"]][[i]][(x1-200):(x1+200),(y1-200):(y1+200)],
            pacientes.semla@tools[["Staffli"]]@rasterlists[["transformed"]][[j]][(x1-200):(x1+200),(y1-200):(y1+200)], 
            listaCoordenadasNEW, c(i,j))
        print(parameters)
        imagescompare <- paste0("comparing_", as.character(i), '_', as.character(j))
        Evaluation$transformed$sameRegion_samePoint[[imagescompare]] <- parameters
        
      }}}}

# Evaluate transformed images (common region, different point)
for ( i in 1:length(listaCoordenadasNEW)) {
  if (i != length(listaCoordenadasNEW)) {
    for ( j in i:length(listaCoordenadasNEW))  {
      if (i != j) {
        x1 <- listaCoordenadasNEW[[i]]$x[[5]]
        y1 <- listaCoordenadasNEW[[i]]$y[[5]]
        x2 <- listaCoordenadasNEW[[j]]$x[[5]]
        y2 <- listaCoordenadasNEW[[j]]$y[[5]]
        print(paste0("Evaluation of the alignment between images ", as.character(i), " and ", as.character(j), " selecting an area from a reference point in each image."))
        
        # Evaluate the alignment and store parameters
        parameters <-
          evalAlign(
            pacientes.semla@tools[["Staffli"]]@rasterlists[["transformed"]][[i]][(x1-200):(x1+200),(y1-200):(y1+200)],
            pacientes.semla@tools[["Staffli"]]@rasterlists[["transformed"]][[j]][(x2-200):(x2+200),(y2-200):(y2+200)],
            listaCoordenadasNEW, c(i,j))
        print(parameters)
        imagescompare <- paste0("comparing_", as.character(i), '_', as.character(j))
        Evaluation$transformed$commonRegion_differentPoint[[imagescompare]] <- parameters
      }}}}


# Update the transformed images in the original Seurat object
pacientes.seurat <- UpdateSeuratFromSemla(pacientes.semla, image_use = "transformed")

############################################################################### Visualize

# Visualize transformed images with spatial dimensions
palign <- SpatialDimPlot(pacientes.seurat, alpha = 0.5, crop = FALSE, ncol = 2,
                         images = c("slice1", "slice2", "slice3", "slice4")) & 
  theme(legend.position = "none",)

# Visualize original images
porig <- SpatialDimPlot(pacientes.seurat, alpha = 0.5, crop = FALSE, ncol = 2,
                        images = c("slice1.1", "slice1.2", "slice1.3", "slice1.4")) & 
  theme(legend.position = "none")

ggsave(paste0(saveDir, "/results/patient19_merge_imageJ_region.png"), plot = palign)
ggsave(paste0(saveDir, "/results/patient19_merge_region.png"), plot = porig)

SpatialDimPlot(pacientes.seurat, alpha = 0.5, crop = FALSE, ncol = 2,
               images = c("slice1.3", "slice3")) & 
  theme(legend.position = "none",)

############################################################################### 

pmt <- SpatialPlot(object = pacientes.seurat, images = c("slice1", "slice2", "slice3", "slice4"), 
                   features = "percent.butterfly", alpha = 1, crop = FALSE, ncol = 2) & 
  theme(legend.position = "top")
ggsave(paste0(saveDir, "/results/patient19_merge_imageJ_genesmt.png"), plot = pmt)

###############################################################################

saveRDS(pacientes.seurat, paste0(saveDir, "/results/patient19_merge_imageJ.rds"))
#pacientes.seurat <- readRDS(paste0(saveDir, "/results/patient19_merge_imageJ.rds"))
#pacientes.semla <- pacientes.seurat
saveRDS(Evaluation, paste0(saveDir, "/results/patient19_EVALUATION_merge_imageJ.rds"))

# Obtain the individual object of the slide
# Splits the Seurat object 'paciente19.merge' into a list of Seurat objects based on the 'name' metadata
pacientes.seurat.split.orig <- SplitObject(paciente19.merge, split.by = "name")

# Generate a spatial dimensional plot for the second object in the split list
# The alpha parameter controls the transparency of the plot, and crop determines whether to crop the image
SpatialDimPlot(pacientes.seurat.split.orig[[2]], alpha = 0.5, crop = FALSE, ncol = 2) & 
  theme(legend.position = "none",)

# Create a copy of the split objects for transformed images
pacientes.seurat.split.trans <- pacientes.seurat.split.orig

# Assign transformed images to the corresponding objects in the new split list
pacientes.seurat.split.trans[[2]]@images$slice1.2 <- pacientes.seurat@images$slice2
pacientes.seurat.split.trans[[3]]@images$slice1.3 <- pacientes.seurat@images$slice3
pacientes.seurat.split.trans[[4]]@images$slice1.4 <- pacientes.seurat@images$slice4

# Generate spatial plots for the "percent.butterfly" feature for the transformed and original objects
p1 <- SpatialPlot(object = pacientes.seurat.split.trans[[3]], 
                  features = "percent.butterfly", alpha = 1, crop = FALSE, ncol = 1) & 
  theme(legend.position = "none")

p2 <- SpatialPlot(object = pacientes.seurat.split.orig[[3]], 
                  features = "percent.butterfly", alpha = 1, crop = FALSE, ncol = 1) & 
  theme(legend.position = "none")

# Print both plots side by side
print(p1 + p2)

# Loop through each of the split original objects starting from the second one
for (i in 2:length(pacientes.seurat.split.orig)) {
  listaObjDeconv <- list()
  
  # Create a list to hold reference and problematic objects for deconvolution
  listaObjDeconv[["reference"]] <- pacientes.seurat.split.orig[[1]]
  listaObjDeconv[["problemNOTalign"]] <- pacientes.seurat.split.orig[[i]]
  listaObjDeconv[["problemYESalign"]] <- pacientes.seurat.split.trans[[i]]
  
  # Update the name in the metadata for each object to reflect its status
  listaObjDeconv[["reference"]]@meta.data$name <- gsub(listaObjDeconv[["reference"]]@meta.data$name[[1]], "reference", listaObjDeconv[["reference"]]@meta.data$name)
  listaObjDeconv[["problemNOTalign"]]@meta.data$name <- gsub(listaObjDeconv[["problemNOTalign"]]@meta.data$name[[1]], "problemNOTalign", listaObjDeconv[["problemNOTalign"]]@meta.data$name)
  listaObjDeconv[["problemYESalign"]]@meta.data$name <- gsub(listaObjDeconv[["problemYESalign"]]@meta.data$name[[1]], "problemYESalign", listaObjDeconv[["problemYESalign"]]@meta.data$name)
  
  savename <- paste0("/results/patient19_merge_imageJ_list_im", i, ".rds")
  saveRDS(listaObjDeconv, paste0(saveDir, savename))
}

############################################################################### Statistical part
Evaluation <- readRDS(paste0(saveDir, "/results/patient19_EVALUATION_merge_imageJ.rds"))

# Initialize vectors to store MSE, SSIM, and Euclidean values for original and transformed images
MSE_controlPos <- c()
MSE_controlNeg <- c()
MSE_controlMov <- c()
MSE_original <- c()
MSE_transformado <- c()
MSE_gray_controlPos <- c()
MSE_gray_controlNeg <- c()
MSE_gray_controlMov <- c()
MSE_gray_original <- c()
MSE_gray_transformado <- c()
SSIM_controlPos <- c()
SSIM_controlNeg <- c()
SSIM_controlMov <- c()
SSIM_original <- c()
SSIM_transformado <- c()
Eucl_original <- c()
Eucl_transformado <- c()

# Loop through the first three images to extract MSE, SSIM, and Euclidean values
for (i in 1:3) {
  MSE_original <- c(MSE_original, Evaluation$original$sameRegion_samePoint[[i]]$mse_value)
  MSE_transformado <- c(MSE_transformado, Evaluation$transformed$sameRegion_samePoint[[i]]$mse_value)
  MSE_gray_original <- c(MSE_gray_original, Evaluation$original$sameRegion_samePoint[[i]]$mse_gray_value)
  MSE_gray_transformado <- c(MSE_gray_transformado, Evaluation$transformed$sameRegion_samePoint[[i]]$mse_gray_value)
  SSIM_original <- c(SSIM_original, Evaluation$original$sameRegion_samePoint[[i]]$ssim_value)
  SSIM_transformado <- c(SSIM_transformado, Evaluation$transformed$sameRegion_samePoint[[i]]$ssim_value)
  Eucl_original <- c(Eucl_original, Evaluation$original$sameRegion_samePoint[[i]]$Eucl_value)
  Eucl_transformado <- c(Eucl_transformado, Evaluation$transformed$sameRegion_samePoint[[i]]$Eucl_value)
}

# Loop through the control data to extract MSE and SSIM values
for (i in 2:4) {
  MSE_controlPos <- c(MSE_controlPos, Evaluation$control[[i]]$`Positive control solution`$mse_value)
  MSE_controlNeg <- c(MSE_controlNeg, Evaluation$control[[i]]$`Negative control solution`$mse_value)
  MSE_controlMov <- c(MSE_controlMov, Evaluation$control[[i]]$`Movement control solution`$mse_value)
  MSE_gray_controlPos <- c(MSE_gray_controlPos, Evaluation$control[[i]]$`Positive control solution`$mse_gray_value)
  MSE_gray_controlNeg <- c(MSE_gray_controlNeg, Evaluation$control[[i]]$`Negative control solution`$mse_gray_value)
  MSE_gray_controlMov <- c(MSE_gray_controlMov, Evaluation$control[[i]]$`Movement control solution`$mse_gray_value)
  SSIM_controlPos <- c(SSIM_controlPos, Evaluation$control[[i]]$`Positive control solution`$ssim_value)
  SSIM_controlNeg <- c(SSIM_controlNeg, Evaluation$control[[i]]$`Negative control solution`$ssim_value)
  SSIM_controlMov <- c(SSIM_controlMov, Evaluation$control[[i]]$`Movement control solution`$ssim_value)
}

# Create a data frame to hold the values from the same region and same point
data_sameRegion_samePoint <- data.frame(
  Image = c("Image2", "Image3", "Image4"),
  MSE_controlPos = MSE_controlPos,
  MSE_controlNeg = MSE_controlNeg,
  MSE_controlMov = MSE_controlMov,
  MSE_original = MSE_original,
  MSE_transformado = MSE_transformado,
  MSE_gray_controlPos = MSE_gray_controlPos,
  MSE_gray_controlNeg = MSE_gray_controlNeg,
  MSE_gray_controlMov = MSE_gray_controlMov,
  MSE_gray_original = MSE_gray_original,
  MSE_gray_transformado = MSE_gray_transformado,
  SSIM_controlPos = SSIM_controlPos,
  SSIM_controlNeg = SSIM_controlNeg,
  SSIM_controlMov = SSIM_controlMov,
  SSIM_original = SSIM_original,
  SSIM_transformado = SSIM_transformado,
  Eucl_original = Eucl_original,
  Eucl_transformado = Eucl_transformado
)

# Set row names for the data frame
rownames(data_sameRegion_samePoint) <- c("Image2", "Image3", "Image4")
data_sameRegion_samePoint

library(ggplot2)
# Create a scatter plot for the differnet parameters
pMSE <- ggplot(data_sameRegion_samePoint, aes(x = Image)) +
  geom_point(aes(y = MSE_controlPos, color = "MSE Positive Control"), size = 5) +
  #  geom_point(aes(y = MSE_controlNeg, color = "MSE Negative Control"), size = 5) +
  geom_point(aes(y = MSE_controlMov, color = "MSE Movement Control"), size = 5) +
  geom_point(aes(y = MSE_original, color = "MSE original"), size = 5) +
  geom_point(aes(y = MSE_transformado, color = "MSE transformed"), size = 5) +
  labs(x = "", y = "MSE value") +
  scale_color_manual(name = "", values = c("MSE Positive Control" = "orange","MSE Negative Control" = "red", "MSE Movement Control" = "blue", "MSE original" = "green", "MSE transformed" = "purple")) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 2))

pMSEgray <- ggplot(data_sameRegion_samePoint, aes(x = Image)) +
  geom_point(aes(y = MSE_gray_controlPos, color = "MSE gray Positive Control"), size = 5) +
  #  geom_point(aes(y = MSE_gray_controlNeg, color = "MSE gray Negative Control"), size = 5) +
  geom_point(aes(y = MSE_gray_controlMov, color = "MSE gray Movement Control"), size = 5) +
  geom_point(aes(y = MSE_gray_original, color = "MSE gray Original"), size = 5) +
  geom_point(aes(y = MSE_gray_transformado, color = "MSE gray Transformed"), size = 5) +
  labs(x = "", y = "MSE gray value") +
  scale_color_manual(name = "", values = c("MSE gray Positive Control" = "orange", "MSE gray Negative Control" = "red", "MSE gray Movement Control" = "blue", "MSE gray Original" = "green", "MSE gray Transformed" = "purple")) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 2))

pSSIM <- ggplot(data_sameRegion_samePoint, aes(x = Image)) +
  geom_point(aes(y = SSIM_controlPos, color = "SSIM Positive Control"), size = 5) +
  geom_point(aes(y = SSIM_controlNeg, color = "SSIM Negative Control"), size = 5) +
  geom_point(aes(y = SSIM_controlMov, color = "SSIM Movement Control"), size = 5) +
  geom_point(aes(y = SSIM_original, color = "SSIM Original"), size = 5) +
  geom_point(aes(y = SSIM_transformado, color = "SSIM Transformed"), size = 5) +
  labs(x = "", y = "SSIM value") +
  scale_color_manual(name = "", values = c("SSIM Positive Control" = "orange", "SSIM Negative Control" = "red", "SSIM Movement Control" = "blue", "SSIM Original" = "green", "SSIM Transformed" = "purple")) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 2))

pEucl <- ggplot(data_sameRegion_samePoint, aes(x = Image)) +
  geom_point(aes(y = Eucl_original, color = "Euclidean Original"), size = 5) +
  geom_point(aes(y = Eucl_transformado, color = "Euclidean Transformed"), size = 5) +
  labs(x = "", y = "Euclidean value") +
  scale_color_manual(name = "", values = c("Euclidean Original" = "green", "Euclidean Transformed" = "purple")) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "bottom"
  )

ggsave(paste0(saveDir, "/results/patient19_merge_imageJ_evalMSE.png"), plot = pMSE)
ggsave(paste0(saveDir, "/results/patient19_merge_imageJ_evalMSEgray.png"), plot = pMSEgray)
ggsave(paste0(saveDir, "/results/patient19_merge_imageJ_evalSSIM.png"), plot = pSSIM)
ggsave(paste0(saveDir, "/results/patient19_merge_imageJ_evalEucl.png"), plot = pEucl)

# Calculate the mean and standard deviation
mean_values <- colMeans(data_sameRegion_samePoint[, -1])  # Ignorar la columna "Image"
sd_values <- apply(data_sameRegion_samePoint[, -1], 2, sd)  # Ignorar la columna "Image"

# Create a new dataframe for statistics
stats_df <- data.frame(
  Parameter = c("MSE Movement Control","MSE Original","MSE Tranformed",
                "MSE gray Movement Control","MSE gray Original","MSE gray Transformed",
                "SSIM Movement Control","SSIM Original","SSIM Transformed",
                "Euclidean Original","Euclidean Transformed"),
  Mean = mean_values,
  SD = sd_values
)

# Perform t-tests
mse_test <- t.test(data_sameRegion_samePoint$MSE_original,
                   data_sameRegion_samePoint$MSE_transformado)
mse_p_value <- mse_test$p.value

ssim_test <- t.test(data_sameRegion_samePoint$SSIM_original,
                    data_sameRegion_samePoint$SSIM_transformado)
ssim_p_value <- ssim_test$p.value

eucl_test <- t.test(data_sameRegion_samePoint$Eucl_original,
                    data_sameRegion_samePoint$Eucl_transformado)
eucl_p_value <- eucl_test$p.value

mse_gray_test <- t.test(data_sameRegion_samePoint$MSE_gray_original,
                        data_sameRegion_samePoint$MSE_gray_transformado)
mse_gray_p_value <- mse_gray_test$p.value

library(tidyr)

# Transform the data to long format
stats_long <- stats_df %>%
  pivot_longer(cols = c(Mean, SD), names_to = "Statistic", values_to = "Value")

# Create bar plots for evaluations
peval1 <- ggplot(stats_df[2:3,], aes(x = Parameter, y = Mean)) +
  geom_bar(stat = "identity", fill = c("green","purple"), width = 0.4) +  # Barras de la media
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +  # Barras de error
  labs(y = "Mean value", x = "", title = "MSE") +
  scale_x_discrete(labels = c("Original", "Transformed")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 28, face = "bold", hjust = 0.5),,
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 28)  # Tamaño del título del eje y
  ) +
  annotate("text", x = 2, y = stats_df[3,]$Mean + 150,
           label = round(mse_p_value, 4), size = 8)

ggsave(paste0(saveDir, "/results/patient19_merge_imageJ_eval1.png"), plot = peval1)

peval2 <- ggplot(stats_df[5:6,], aes(x = Parameter, y = Mean)) +
  geom_bar(stat = "identity", fill = c("green","purple"), width = 0.4) +  # Barras de la media
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +  # Barras de error
  labs(y = "Mean value", x = "", title = "MSE gray") +
  scale_x_discrete(labels = c("Original", "Transformed")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 28, face = "bold", hjust = 0.5),,
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 28)  # Tamaño del título del eje y
  ) +
  annotate("text", x = 2, y = stats_df[6,]$Mean + 150,
           label = round(mse_gray_p_value, 4), size = 8)

ggsave(paste0(saveDir, "/results/patient19_merge_imageJ_eval2.png"), plot = peval2)

peval3 <- ggplot(stats_df[8:9,], aes(x = Parameter, y = Mean)) +
  geom_bar(stat = "identity", fill = c("green","purple"), width = 0.4) +  # Barras de la media
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +  # Barras de error
  labs(y = "Mean value", x = "", title = "SSIM") +
  scale_x_discrete(labels = c("Original", "Transformed")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 28, face = "bold", hjust = 0.5),,
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 28)  # Tamaño del título del eje y
  ) +
  annotate("text", x = 2, y = stats_df[9,]$Mean + 0.05,
           label = round(ssim_p_value, 4), size = 8)

ggsave(paste0(saveDir, "/results/patient19_merge_imageJ_eval3.png"), plot = peval3)

peval4 <- ggplot(stats_df[10:11,], aes(x = Parameter, y = Mean)) +
  geom_bar(stat = "identity", fill = c("green","purple"), width = 0.4) +  # Barras de la media
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +  # Barras de error
  labs(y = "Mean value", x = "", title = "Euclidean distance") +
  scale_x_discrete(labels = c("Original", "Transformed")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 28, face = "bold", hjust = 0.5),,
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 28)  # Tamaño del título del eje y
  ) +
  annotate("text", x = 2, y = stats_df[11,]$Mean + 150,
           label = round(eucl_p_value, 4), size = 8)

ggsave(paste0(saveDir, "/results/patient19_merge_imageJ_eval4.png"), plot = peval4)
