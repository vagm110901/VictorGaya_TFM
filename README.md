# Development of image alignment methods for Spatial Transcriptomics data of Spinal Cord

This repository contains the scripts used for the Masther's Thesis titled _"Development of image alignment methods for Spatial Transcriptomics data of Spinal Cord"_ done by **Víctor Álvaro Gaya Martín**, supervissed by the external tutor **Ana Victoria Conesa Cegarra**, having as UV tutor Juli Peretó Magraner, for the Masther in Bioinformatics - UV.

  * Inside the _dataPreparation_ directory are the scripts used to read spatial transcriptomics and single nucleus data from raw files to prepare the Seurat objects for the analysis.
  * Inside the _mirrorTransformation_ directory are the scripts used to do the probes for the specular transformation.
  * Inside the _alignmentMethods_ directory are the directories for each different approach (_Geomtric_, _Procrustes_ and _imageJ_) and a directory for evaluation scripts.
    * Inside the _Geometric_ directory are the scripts which use the **Geometric Transformation Estimation Model** for the alignment.
    * Inside the _Procrustes_ directory are the scripts which use **IMIFA::Procrustes function** for the alignment.
    * Inside the _imageJ_ directory are the scripts which use the **ImageJ plugin Register Virtual Stack Slices** for the alignment.
    * Inside the _evaluation_ directory are the scripts for the statitical and biological evaluation of the alignment.
