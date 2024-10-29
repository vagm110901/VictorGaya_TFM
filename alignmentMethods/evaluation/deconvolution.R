# Load required libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)

# Define directories for data and output
carpetaData <- "/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/02.Seurat/Yadav2023_adult_human_spinal_cord/snRNAseq_results"
saveDir <- "/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/04.ImageAlingment"
saveDir <- "/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/05.ImageJ"

# Load pre-processed Spatial Transcriptomics (ST) and snRNAseq data
listaObjST <- readRDS(paste0(saveDir, "/results/patient19_merge_alignTR_list_im2.rds"))
listaObjST <- readRDS(paste0(saveDir, "/results/patientMIX_merge_alignTR_list_im2.rds"))
singleNucleus <- readRDS(paste0(carpetaData, "/snRNA_integratedCCA.rds"))

# Define subset range for spatial coordinates in the ST object
objectST <- listaObjST[[1]]
imagerowmin <- min(objectST@images[[1]]@coordinates$imagerow[objectST@images[[1]]@coordinates$row > 20 & 
                                                               objectST@images[[1]]@coordinates$row < 43])
imagerowmax <- max(objectST@images[[1]]@coordinates$imagerow[objectST@images[[1]]@coordinates$row > 20 & 
                                                               objectST@images[[1]]@coordinates$row < 43])
imagecolmin <- min(objectST@images[[1]]@coordinates$imagecol[objectST@images[[1]]@coordinates$col > 20 & 
                                                               objectST@images[[1]]@coordinates$col < 61])
imagecolmax <- max(objectST@images[[1]]@coordinates$imagecol[objectST@images[[1]]@coordinates$col > 20 & 
                                                               objectST@images[[1]]@coordinates$col < 61])

# Initialize lists for zoomed ST objects and images
listaObjZoom <- list()
listaImages <- list()

# Loop through ST objects and create subsetted, zoomed objects
for (i in seq_along(listaObjST)) {
  objectST <- listaObjST[[i]]
  
  # Define cells within subset area
  col <- objectST@images[[1]]@coordinates$imagecol
  row <- objectST@images[[1]]@coordinates$imagerow
  grupo <- rownames(objectST@images[[1]]@coordinates[which(col >= imagecolmin & col <= imagecolmax & 
                                                             row >= imagerowmin & row <= imagerowmax),])
  
  # Subset and re-normalize each ST object
  objectST.zoom <- subset(objectST, cells = grupo, invert = FALSE) %>%
    SCTransform(assay = "Spatial", verbose = FALSE) %>%
    RunPCA(verbose = FALSE)
  objname <- names(listaObjST[i])
  listaObjZoom[[objname]] <- objectST.zoom
}

# Plot the zoomed SpatialDimPlots for each ST object
m1 <- SpatialDimPlot(listaObjZoom[[1]], alpha = 0.5, crop = FALSE) & theme(legend.position = "none")
m2 <- SpatialDimPlot(listaObjZoom[[2]], alpha = 0.5, crop = FALSE) & theme(legend.position = "none")
m3 <- SpatialDimPlot(listaObjZoom[[3]], alpha = 0.5, crop = FALSE) & theme(legend.position = "none")
print(m1 + m2 + m3)

# Transfer cell type annotations from snRNAseq to the zoomed ST objects
DimPlot(singleNucleus, group.by = "seurat_clusters", label = TRUE)
DimPlot(singleNucleus, group.by = "cell.types", label = TRUE)

listaObjAnnot <- list()
for (i in seq_along(listaObjST)) {
  objectST <- listaObjZoom[[i]]
  anchors <- FindTransferAnchors(reference = singleNucleus, query = objectST, normalization.method = "SCT")
  predictions.assay <- TransferData(anchorset = anchors, refdata = singleNucleus$cell.types, 
                                    prediction.assay = TRUE, weight.reduction = objectST[["pca"]], 
                                    dim = 1:30, k.weight = 30)
  objectST[["predictions"]] <- predictions.assay
  DefaultAssay(objectST) <- "predictions"
  
  objname <- names(listaObjST[i])
  listaObjAnnot[[objname]] <- objectST
}

# Save annotated ST objects
saveRDS(listaObjAnnot, paste0(saveDir, "/patient19_merge_alignTR_list_im2_annot.rds"))

# Load annotated ST objects for visualization
listaObjAnnot <- readRDS(paste0(saveDir, "/results/patientMIX_merge_imageJ_list_im2_annot.rds"))

# Generate SpatialFeaturePlots for cell types in each annotated ST object
q1 <- SpatialFeaturePlot(listaObjAnnot[[1]], features = "neurons", pt.size.factor = 1, ncol = 2, crop = FALSE) + 
  scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"), 
                       limits = c(0, 1), 
                       oob = scales::squish)
q2 <- SpatialFeaturePlot(listaObjAnnot[[1]], features = "oligodendrocytes", pt.size.factor = 1, ncol = 1, crop = FALSE) + 
  scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"), 
                       limits = c(0, 1), 
                       oob = scales::squish)

q3 <- SpatialFeaturePlot(listaObjAnnot[[2]], features = "neurons", pt.size.factor = 1, ncol = 2, crop = FALSE) + 
  scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"), 
                       limits = c(0, 1), 
                       oob = scales::squish)
q4 <- SpatialFeaturePlot(listaObjAnnot[[2]], features = "oligodendrocytes", pt.size.factor = 1, ncol = 1, crop = FALSE) + 
  scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"), 
                       limits = c(0, 1), 
                       oob = scales::squish)

q5 <- SpatialFeaturePlot(listaObjAnnot[[3]], features = "neurons", pt.size.factor = 1, ncol = 2, crop = FALSE) + 
  scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"), 
                       limits = c(0, 1), 
                       oob = scales::squish)
q6 <- SpatialFeaturePlot(listaObjAnnot[[3]], features = "oligodendrocytes", pt.size.factor = 1, ncol = 1, crop = FALSE) + 
  scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"), 
                       limits = c(0, 1), 
                       oob = scales::squish)

# Save plots to PNG
png(paste0(saveDir, "/results/patientMIX_merge_imageJ_deconvORIG.png"))
q1 + q2
dev.off()

png(paste0(saveDir, "/results/patientMIX_merge_imageJ_deconv2NO.png"))
q3 + q4
dev.off()

png(paste0(saveDir, "/results/patientMIX_merge_imageJ_deconv2SI.png"))
q5 + q6
dev.off()

# Deconvolution analysis for all cell types for the non-aligned object
cell.types <- c("astrocytes", "microglia", "OPCs", "endothelial cells", "neurons", 
                "pericytes+meningeal+Swchann cells", "lymphocytes", "oligodendrocytes")

listaGraphs <- list()
for (type in cell.types) {
  listaGraphs[[type]] <- SpatialFeaturePlot(listaObjAnnot[[1]], features = type, pt.size.factor = 1, crop = FALSE) +
    scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"), limits = c(0, 1), oob = scales::squish)
}

# Save deconvolution plots for various cell types to PNG files
png(paste0(saveDir, "/results/patient19_merge_DECONV_original1.png"))
( listaGraphs[[1]] + listaGraphs[[2]] ) 
dev.off()
png(paste0(saveDir, "/results/patient19_merge_DECONV_original2.png"))
( listaGraphs[[3]] + listaGraphs[[4]] ) 
dev.off()
png(paste0(saveDir, "/results/patient19_merge_DECONV_original3.png"))
( listaGraphs[[5]] + listaGraphs[[6]] ) 
dev.off()
png(paste0(saveDir, "/results/patient19_merge_DECONV_original4.png"))
( listaGraphs[[7]] + listaGraphs[[8]] )
dev.off()
