####### Integration of all datasets for Simpson et al., 2024 #######

### Seurat script for single cell analysis 
{
  # Load libraries required for the analysis
  library(glmGamPoi)  # For GLM-based normalization of gene expression data
  library(Seurat)     # Core library for single cell genomics analysis
  library(SeuratDisk) # For interaction with Seurat objects on disk
  library(patchwork)  # For arranging ggplot2 plots
  library(dplyr)      # For data manipulation
  library(cowplot)    # For additional plotting features
  library(ggplot2)    # For data visualization
  library(leidenbase) # For clustering algorithms
  library(SeuratWrappers) # Extensions for the Seurat package
  library(Matrix)     # For sparse and other matrix methods
  library(tidyverse)  # For data manipulation and visualization
  library(DoubletFinder) # For doublet identification and removal
}

# The code below was used for the individual processing of samples before integration, 
# however, given that samples were added and processed sequentially this section contains a large amount 
# of repetitive code consider merging objects and processing as a list

# Setwd
setwd("/Users/lukesimpson/Desktop/Integrated analysis")
# Read in data
PigE11.5A.EPS <- Read10X(data.dir = "/Users/lukesimpson/Desktop/New SC RNA seq data 07112022/ds973_1/outs/filtered_feature_bc_matrix")
PigE11.5B.EPS <- Read10X(data.dir = "/Users/lukesimpson/Desktop/New SC RNA seq data 07112022/ds973_2/outs/filtered_feature_bc_matrix")
PigE12A.PS <- Read10X(data.dir = "/Users/lukesimpson/Desktop/New SC RNA seq data 07112022/ds973_3/outs/filtered_feature_bc_matrix")
PigE12B.PS <- Read10X(data.dir = "/Users/lukesimpson/Desktop/All Pig samples/E12 B (DS973_4)/outs/filtered_feature_bc_matrix")
PigE12.5A.PS2 <- Read10X(data.dir = "/Users/lukesimpson/Desktop/New SC RNA seq data 07112022/ds973_18/outs/filtered_feature_bc_matrix")
PigE12.5B.PS2 <- Read10X(data.dir = "/Users/lukesimpson/Desktop/New SC RNA seq data 07112022/ds973_19/outs/filtered_feature_bc_matrix")
PigE13A.NP <- Read10X(data.dir = "/Users/lukesimpson/Desktop/All Pig samples/E13.5 A (DS973_8)/outs/filtered_feature_bc_matrix")
PigE13B.NP <- Read10X(data.dir = "/Users/lukesimpson/Desktop/New SC RNA seq data 07112022/DS1149_1  Filtered_feature_bc_matrix")
PigE13C.NP <- Read10X(data.dir = "/Users/lukesimpson/Desktop/New SC RNA seq data 07112022/DS1149_2 filtered_feature_bc_matrix")
PigE13D.NP <- Read10X(data.dir = "/Users/lukesimpson/Desktop/New SC RNA seq data 07112022/DS1149_5 filtered_feature_bc_matrix")
PigE13E.NP <- Read10X(data.dir = "/Users/lukesimpson/Desktop/New SC RNA seq data 07112022/DS1149_6/outs/filtered_feature_bc_matrix")
PigE13.5A.NG <- Read10X(data.dir = "/Users/lukesimpson/Desktop/All Pig samples/Neural groove A (DS973_20)/outs/filtered_feature_bc_matrix")
PigE13.5B.NG <- Read10X(data.dir = "/Users/lukesimpson/Desktop/All Pig samples/Neural groove B (DS973_21)/outs/filtered_feature_bc_matrix")
PigE13.5C.NG <- Read10X(data.dir = "/Users/lukesimpson/Desktop/New SC RNA seq data 07112022/DS1149_7 filtered_feature_bc_matrix")
PigE13.5D.NG <- Read10X(data.dir = "/Users/lukesimpson/Desktop/New SC RNA seq data 07112022/DS1149_8 filtered_feature_bc_matrix")
PigE14A.6SOM <- Read10X(data.dir = "/Users/lukesimpson/Desktop/New SC RNA seq data 07112022/ds973_9/outs/filtered_feature_bc_matrix")
PigE14B.6SOM <- Read10X(data.dir = "/Users/lukesimpson/Desktop/New SC RNA seq data 07112022/ds973_10/outs/filtered_feature_bc_matrix")
PigE14.5A.8SOM <- Read10X(data.dir = "/Users/lukesimpson/Desktop/New SC RNA seq data 07112022/ds973_15/outs/filtered_feature_bc_matrix")
PigE14.5B.8SOM <- Read10X(data.dir = "/Users/lukesimpson/Desktop/All Pig samples/E14 C (DS973_5)/outs/filtered_feature_bc_matrix")
PigE14.5C.8SOM <- Read10X(data.dir = "/Users/lukesimpson/Desktop/All Pig samples/E14 Posterior only (DS973_6)/outs/filtered_feature_bc_matrix")
PigE14.5D.8SOM <- Read10X(data.dir = "/Users/lukesimpson/Desktop/All Pig samples/E14 Posterior only (DS973_7)/outs/filtered_feature_bc_matrix")
PigE15A.10SOM <- Read10X(data.dir = "/Users/lukesimpson/Desktop/All Pig samples/E15 (DS973_11)/outs/filtered_feature_bc_matrix")
PigE15B.10SOM <- Read10X(data.dir = "/Users/lukesimpson/Desktop/New SC RNA seq data 07112022/DS1149_3 filtered_feature_bc_matrix")
PigE15C.10SOM <- Read10X(data.dir = "/Users/lukesimpson/Desktop/New SC RNA seq data 07112022/ds1149_4/outs/filtered_feature_bc_matrix")

## Create the Seurat objects, min.features was chosen based on the distributions of the samples (NOTE usually they have a binominal distribution and low quality cells were below 1750 in most cases)
PigE11.5A.EPS <- CreateSeuratObject(counts = PigE11.5A.EPS, project ="E11.5A.EPS" , min.cells = 1, min.features = 1750, assay = "RNA")
PigE11.5B.EPS <- CreateSeuratObject(counts = PigE11.5B.EPS, project ="E11.5B.EPS" , min.cells = 1, min.features = 1750, assay = "RNA")
PigE12A.PS <- CreateSeuratObject(counts = PigE12A.PS,project ="E12A.PS" , min.cells = 1, min.features = 1750, assay = "RNA")
PigE12B.PS <- CreateSeuratObject(counts = PigE12B.PS,project ="E12B.PS" , min.cells = 1, min.features = 1750, assay = "RNA")
PigE12.5A.PS2 <- CreateSeuratObject(counts = PigE12.5A.PS2, project ="E12.5A.PS2" ,min.cells = 1, min.features = 1750, assay = "RNA")
PigE12.5B.PS2 <- CreateSeuratObject(counts = PigE12.5B.PS2,project ="E12.5B.PS2" , min.cells = 1, min.features = 1750, assay = "RNA")
PigE13A.NP <- CreateSeuratObject(counts = PigE13A.NP, project ="E13A.NP" ,min.cells = 1, min.features = 1750, assay = "RNA")
PigE13B.NP <- CreateSeuratObject(counts = PigE13B.NP,project ="E13B.NP" , min.cells = 1, min.features = 1750, assay = "RNA")
PigE13C.NP <- CreateSeuratObject(counts = PigE13C.NP,project ="E13C.NP" , min.cells = 1, min.features = 1750, assay = "RNA")
PigE13D.NP <- CreateSeuratObject(counts = PigE13D.NP,project ="E13D.NP" , min.cells = 1, min.features = 1750, assay = "RNA")
PigE13E.NP <- CreateSeuratObject(counts = PigE13E.NP,project ="E13E.NP" , min.cells = 1, min.features = 1750, assay = "RNA")
PigE13.5A.NG <- CreateSeuratObject(counts = PigE13.5A.NG,project ="E13.5A.NG" , min.cells = 1, min.features = 1750, assay = "RNA")
PigE13.5B.NG <- CreateSeuratObject(counts = PigE13.5B.NG,project ="E13.5B.NG" , min.cells = 1, min.features = 1750, assay = "RNA")
PigE13.5C.NG <- CreateSeuratObject(counts = PigE13.5C.NG,project ="E13.5C.NG" , min.cells = 1, min.features = 1750, assay = "RNA")
PigE13.5D.NG <- CreateSeuratObject(counts = PigE13.5D.NG,project ="E13.5D.NG" , min.cells = 1, min.features = 1750, assay = "RNA")
PigE14A.6SOM <- CreateSeuratObject(counts = PigE14A.6SOM,project ="E14A.6SOM" , min.cells = 1, min.features = 1750, assay = "RNA")
PigE14B.6SOM <- CreateSeuratObject(counts = PigE14B.6SOM, project ="E14B.6SOM" ,min.cells = 1, min.features = 1750, assay = "RNA")
PigE14.5A.8SOM <- CreateSeuratObject(counts = PigE14.5A.8SOM,project ="E14.5A.8SOM" , min.cells = 1, min.features = 1750, assay = "RNA")
PigE14.5B.8SOM <- CreateSeuratObject(counts = PigE14.5B.8SOM, project ="E14.5B.8SOM" ,min.cells = 1, min.features = 1750, assay = "RNA")
PigE14.5C.8SOM <- CreateSeuratObject(counts = PigE14.5C.8SOM, project ="E14.5C.8SOM" ,min.cells = 1, min.features = 1750, assay = "RNA")
PigE14.5D.8SOM <- CreateSeuratObject(counts = PigE14.5D.8SOM,project ="E14.5D.8SOM" ,min.cells = 1, min.features = 1750, assay = "RNA")
PigE15A.10SOM <- CreateSeuratObject(counts = PigE15A.10SOM, project ="E15A.10SOM" ,min.cells = 1, min.features = 1750, assay = "RNA")
PigE15B.10SOM <- CreateSeuratObject(counts = PigE15B.10SOM,project ="E15B.10SOM" , min.cells = 1, min.features = 1750, assay = "RNA")
PigE15C.10SOM <- CreateSeuratObject(counts = PigE15C.10SOM,project ="E15C.10SOM" , min.cells = 1, min.features = 1750, assay = "RNA")

# Identify apoptotic cells by % Mitochondrial gene expression (NOTE that I do not use this to directly fiter out cells)
Mitochondrial_genes_pig_orthos <- c("ATP6",
                                    "ATP8",
                                    "COX1",
                                    "COX2",
                                    "COX3",
                                    "CYTB",
                                    "ND2",
                                    "ND3",
                                    "ND4",
                                    "ND4L",
                                    "ND5",
                                    "ND6")

Y_chromosome_genes_pig_orthos <- c("SRY",
                                   "ZFY",
                                   "RPS4Y1",
                                   "AMELY",
                                   "TBL1Y",
                                   "PCDH11Y",
                                   "TGIF2LY",
                                   "TSPY1",
                                   "USP9Y",
                                   "DDX3Y",
                                   "UTY",
                                   "TB4Y",
                                   "RPS4Y2",
                                   "EIF1AY",
                                   "KDM5D",
                                   "XKRY",
                                   "HSFY1",
                                   "PRY",
                                   "RBMY1A1",
                                   "AZFC",
                                   "DAZ1",
                                   "CDY1",
                                   "VCY1",
                                   "VCY2",
                                   "DAZ2",
                                   "DAZ2",
                                   "TSPY2",
                                   "DAZ4",
                                   "CDY2",
                                   "PRY2",
                                   "HSFY2")

#Calculate percentage MT genes
{
  PigE11.5A.EPS[["percent.RP"]] <- PercentageFeatureSet(PigE11.5A.EPS,pattern = "^RP")
  PigE11.5B.EPS[["percent.RP"]] <- PercentageFeatureSet(PigE11.5B.EPS,pattern = "^RP")
  PigE12A.PS[["percent.RP"]] <- PercentageFeatureSet(PigE12A.PS,pattern = "^RP")
  PigE12B.PS[["percent.RP"]] <- PercentageFeatureSet(PigE12B.PS,pattern = "^RP")
  PigE12.5A.PS2[["percent.RP"]] <- PercentageFeatureSet(PigE12.5A.PS2, pattern = "^RP")
  PigE12.5B.PS2[["percent.RP"]] <- PercentageFeatureSet(PigE12.5B.PS2,pattern = "^RP")
  PigE13A.NP[["percent.RP"]] <- PercentageFeatureSet(PigE13A.NP,pattern = "^RP")
  PigE13B.NP[["percent.RP"]] <- PercentageFeatureSet(PigE13B.NP, pattern = "^RP")
  PigE13C.NP[["percent.RP"]] <- PercentageFeatureSet(PigE13C.NP,pattern = "^RP")
  PigE13D.NP[["percent.RP"]] <- PercentageFeatureSet(PigE13D.NP,pattern = "^RP")
  PigE13E.NP[["percent.RP"]] <- PercentageFeatureSet(PigE13E.NP, pattern = "^RP")
  PigE13.5A.NG[["percent.RP"]] <- PercentageFeatureSet(PigE13.5A.NG, pattern = "^RP")
  PigE13.5B.NG[["percent.RP"]] <- PercentageFeatureSet(PigE13.5B.NG, pattern = "^RP")
  PigE13.5C.NG[["percent.RP"]] <- PercentageFeatureSet(PigE13.5C.NG, pattern = "^RP")
  PigE13.5D.NG[["percent.RP"]] <- PercentageFeatureSet(PigE13.5D.NG, pattern = "^RP")
  PigE14A.6SOM[["percent.RP"]] <- PercentageFeatureSet(PigE14A.6SOM, pattern = "^RP")
  PigE14B.6SOM[["percent.RP"]] <- PercentageFeatureSet(PigE14B.6SOM, pattern = "^RP")
  PigE14.5A.8SOM[["percent.RP"]] <- PercentageFeatureSet(PigE14.5A.8SOM,pattern = "^RP")
  PigE14.5B.8SOM[["percent.RP"]] <- PercentageFeatureSet(PigE14.5B.8SOM,pattern = "^RP")
  PigE14.5C.8SOM[["percent.RP"]] <- PercentageFeatureSet(PigE14.5C.8SOM,pattern = "^RP")
  PigE14.5D.8SOM[["percent.RP"]] <- PercentageFeatureSet(PigE14.5D.8SOM, pattern = "^RP")
  PigE15A.10SOM[["percent.RP"]] <- PercentageFeatureSet(PigE15A.10SOM,pattern = "^RP")
  PigE15B.10SOM[["percent.RP"]] <- PercentageFeatureSet(PigE15B.10SOM,pattern = "^RP")
  PigE15C.10SOM[["percent.RP"]] <- PercentageFeatureSet(PigE15C.10SOM,pattern = "^RP")
}

{
  PigE11.5A.EPS[["percent.Ychromo"]] <- PercentageFeatureSet(PigE11.5A.EPS, features = Y_chromosome_genes_pig_orthos)
  PigE11.5B.EPS[["percent.Ychromo"]] <- PercentageFeatureSet(PigE11.5B.EPS,features = Y_chromosome_genes_pig_orthos)
  PigE12A.PS[["percent.Ychromo"]] <- PercentageFeatureSet(PigE12A.PS,features = Y_chromosome_genes_pig_orthos)
  PigE12B.PS[["percent.Ychromo"]] <- PercentageFeatureSet(PigE12B.PS,features = Y_chromosome_genes_pig_orthos)
  PigE12.5A.PS2[["percent.Ychromo"]] <- PercentageFeatureSet(PigE12.5A.PS2, features =Y_chromosome_genes_pig_orthos)
  PigE12.5B.PS2[["percent.Ychromo"]] <- PercentageFeatureSet(PigE12.5B.PS2,features = Y_chromosome_genes_pig_orthos)
  PigE13A.NP[["percent.Ychromo"]] <- PercentageFeatureSet(PigE13A.NP,features = Y_chromosome_genes_pig_orthos)
  PigE13B.NP[["percent.Ychromo"]] <- PercentageFeatureSet(PigE13B.NP, features =Y_chromosome_genes_pig_orthos)
  PigE13C.NP[["percent.Ychromo"]] <- PercentageFeatureSet(PigE13C.NP,features = Y_chromosome_genes_pig_orthos)
  PigE13D.NP[["percent.Ychromo"]] <- PercentageFeatureSet(PigE13D.NP,features = Y_chromosome_genes_pig_orthos)
  PigE13E.NP[["percent.Ychromo"]] <- PercentageFeatureSet(PigE13E.NP, features =Y_chromosome_genes_pig_orthos)
  PigE13.5A.NG[["percent.Ychromo"]] <- PercentageFeatureSet(PigE13.5A.NG, features =Y_chromosome_genes_pig_orthos)
  PigE13.5B.NG[["percent.Ychromo"]] <- PercentageFeatureSet(PigE13.5B.NG, features =Y_chromosome_genes_pig_orthos)
  PigE13.5C.NG[["percent.Ychromo"]] <- PercentageFeatureSet(PigE13.5C.NG, features =Y_chromosome_genes_pig_orthos)
  PigE13.5D.NG[["percent.Ychromo"]] <- PercentageFeatureSet(PigE13.5D.NG, features =Y_chromosome_genes_pig_orthos)
  PigE14A.6SOM[["percent.Ychromo"]] <- PercentageFeatureSet(PigE14A.6SOM, features =Y_chromosome_genes_pig_orthos)
  PigE14B.6SOM[["percent.Ychromo"]] <- PercentageFeatureSet(PigE14B.6SOM, features =Y_chromosome_genes_pig_orthos)
  PigE14.5A.8SOM[["percent.Ychromo"]] <- PercentageFeatureSet(PigE14.5A.8SOM,features = Y_chromosome_genes_pig_orthos)
  PigE14.5B.8SOM[["percent.Ychromo"]] <- PercentageFeatureSet(PigE14.5B.8SOM,features = Y_chromosome_genes_pig_orthos)
  PigE14.5C.8SOM[["percent.Ychromo"]] <- PercentageFeatureSet(PigE14.5C.8SOM,features = Y_chromosome_genes_pig_orthos)
  PigE14.5D.8SOM[["percent.Ychromo"]] <- PercentageFeatureSet(PigE14.5D.8SOM, features =Y_chromosome_genes_pig_orthos)
  PigE15A.10SOM[["percent.Ychromo"]] <- PercentageFeatureSet(PigE15A.10SOM,features = Y_chromosome_genes_pig_orthos)
  PigE15B.10SOM[["percent.Ychromo"]] <- PercentageFeatureSet(PigE15B.10SOM,features = Y_chromosome_genes_pig_orthos)
  PigE15C.10SOM[["percent.Ychromo"]] <- PercentageFeatureSet(PigE15C.10SOM,features = Y_chromosome_genes_pig_orthos)
}

{
  PigE11.5A.EPS[["percent.mt"]] <- PercentageFeatureSet(PigE11.5A.EPS,features = Mitochondrial_genes_pig_orthos)
  PigE11.5B.EPS[["percent.mt"]] <- PercentageFeatureSet(PigE11.5B.EPS, features =Mitochondrial_genes_pig_orthos)
  PigE12A.PS[["percent.mt"]] <- PercentageFeatureSet(PigE12A.PS, features =Mitochondrial_genes_pig_orthos)
  PigE12B.PS[["percent.mt"]] <- PercentageFeatureSet(PigE12B.PS, features =Mitochondrial_genes_pig_orthos)
  PigE12.5A.PS2[["percent.mt"]] <- PercentageFeatureSet(PigE12.5A.PS2, features =Mitochondrial_genes_pig_orthos)
  PigE12.5B.PS2[["percent.mt"]] <- PercentageFeatureSet(PigE12.5B.PS2,features = Mitochondrial_genes_pig_orthos)
  PigE13A.NP[["percent.mt"]] <- PercentageFeatureSet(PigE13A.NP,features = Mitochondrial_genes_pig_orthos)
  PigE13B.NP[["percent.mt"]] <- PercentageFeatureSet(PigE13B.NP,features = Mitochondrial_genes_pig_orthos)
  PigE13C.NP[["percent.mt"]] <- PercentageFeatureSet(PigE13C.NP,features = Mitochondrial_genes_pig_orthos)
  PigE13D.NP[["percent.mt"]] <- PercentageFeatureSet(PigE13D.NP,features = Mitochondrial_genes_pig_orthos)
  PigE13E.NP[["percent.mt"]] <- PercentageFeatureSet(PigE13E.NP, features =Mitochondrial_genes_pig_orthos)
  PigE13.5A.NG[["percent.mt"]] <- PercentageFeatureSet(PigE13.5A.NG, features =Mitochondrial_genes_pig_orthos)
  PigE13.5B.NG[["percent.mt"]] <- PercentageFeatureSet(PigE13.5B.NG,features = Mitochondrial_genes_pig_orthos)
  PigE13.5C.NG[["percent.mt"]] <- PercentageFeatureSet(PigE13.5C.NG,features = Mitochondrial_genes_pig_orthos)
  PigE13.5D.NG[["percent.mt"]] <- PercentageFeatureSet(PigE13.5D.NG, features =Mitochondrial_genes_pig_orthos)
  PigE14A.6SOM[["percent.mt"]] <- PercentageFeatureSet(PigE14A.6SOM, features =Mitochondrial_genes_pig_orthos)
  PigE14B.6SOM[["percent.mt"]] <- PercentageFeatureSet(PigE14B.6SOM,features = Mitochondrial_genes_pig_orthos)
  PigE14.5A.8SOM[["percent.mt"]] <- PercentageFeatureSet(PigE14.5A.8SOM,features = Mitochondrial_genes_pig_orthos)
  PigE14.5B.8SOM[["percent.mt"]] <- PercentageFeatureSet(PigE14.5B.8SOM,features = Mitochondrial_genes_pig_orthos)
  PigE14.5C.8SOM[["percent.mt"]] <- PercentageFeatureSet(PigE14.5C.8SOM,features = Mitochondrial_genes_pig_orthos)
  PigE14.5D.8SOM[["percent.mt"]] <- PercentageFeatureSet(PigE14.5D.8SOM,features = Mitochondrial_genes_pig_orthos)
  PigE15A.10SOM[["percent.mt"]] <- PercentageFeatureSet(PigE15A.10SOM,features = Mitochondrial_genes_pig_orthos)
  PigE15B.10SOM[["percent.mt"]] <- PercentageFeatureSet(PigE15B.10SOM,features = Mitochondrial_genes_pig_orthos)
  PigE15C.10SOM[["percent.mt"]] <- PercentageFeatureSet(PigE15C.10SOM,features = Mitochondrial_genes_pig_orthos)
}

### QUALITY CONTROL PLOTS ###
#Now that we can visualise the distibutions we can remove cells that have a very high feature count (potential doublets)
# or remove apoptotic cells which have high mitochondrial gene expression
# PigE13B <- subset(PigE13B, subset = nFeature_RNA >1000 & nFeature_RNA < 7500 & percent.mt < 5)
#E11.5A

VlnPlot(PigE11.5A.EPS, features = "nFeature_RNA", pt.size = 0.1) + NoLegend() #Bimodal distribution the majority of cells have around 800 or 4000 genes detected
VlnPlot(PigE11.5B.EPS, features = "nFeature_RNA", pt.size = 0.1) + NoLegend() #Bimodal distribution the majority of cells have around 800 or 4000 genes detected
VlnPlot(PigE12A.PS, features = "nFeature_RNA", pt.size = 0.1) + NoLegend() #Bimodal distribution the majority of cells have around 800 or 4000 genes detected
VlnPlot(PigE12B.PS, features = "nFeature_RNA", pt.size = 0.1) + NoLegend() #Bimodal distribution the majority of cells have around 800 or 4000 genes detected
VlnPlot(PigE12.5A.PS2, features = "nFeature_RNA", pt.size = 0.1) + NoLegend() #Bimodal distribution the majority of cells have around 800 or 4000 genes detected
VlnPlot(PigE12.5B.PS2, features = "nFeature_RNA", pt.size = 0.1) + NoLegend() #Bimodal distribution the majority of cells have around 800 or 4000 genes detected
VlnPlot(PigE13A.NP, features = "nFeature_RNA", pt.size = 0.1) + NoLegend() #Bimodal distribution the majority of cells have around 800 or 4000 genes detected
VlnPlot(PigE13B.NP, features = "nFeature_RNA", pt.size = 0.1) + NoLegend() #Bimodal distribution the majority of cells have around 800 or 4000 genes detected
VlnPlot(PigE13C.NP, features = "nFeature_RNA", pt.size = 0.1) + NoLegend() #Bimodal distribution the majority of cells have around 800 or 4000 genes detected
VlnPlot(PigE13E.NP, features = "nFeature_RNA", pt.size = 0.1) + NoLegend() #Bimodal distribution the majority of cells have around 800 or 4000 genes detected
VlnPlot(PigE13.5A.NG, features = "nFeature_RNA", pt.size = 0.1) + NoLegend() #Bimodal distribution the majority of cells have around 800 or 4000 genes detected
VlnPlot(PigE13.5B.NG, features = "nFeature_RNA", pt.size = 0.1) + NoLegend() #Bimodal distribution the majority of cells have around 800 or 4000 genes detected
VlnPlot(PigE13.5C.NG, features = "nFeature_RNA", pt.size = 0.1) + NoLegend() #Bimodal distribution the majority of cells have around 800 or 4000 genes detected
VlnPlot(PigE13.5D.NG, features = "nFeature_RNA", pt.size = 0.1) + NoLegend() #Bimodal distribution the majority of cells have around 800 or 4000 genes detected
VlnPlot(PigE14A.6SOM, features = "nFeature_RNA", pt.size = 0.1) + NoLegend() #Bimodal distribution the majority of cells have around 800 or 4000 genes detected
VlnPlot(PigE14B.6SOM, features = "nFeature_RNA", pt.size = 0.1) + NoLegend() #Bimodal distribution the majority of cells have around 800 or 4000 genes detected
VlnPlot(PigE14.5A.8SOM, features = "nFeature_RNA", pt.size = 0.1) + NoLegend() #Bimodal distribution the majority of cells have around 800 or 4000 genes detected
VlnPlot(PigE14.5B.8SOM, features = "nFeature_RNA", pt.size = 0.1) + NoLegend() #Bimodal distribution the majority of cells have around 800 or 4000 genes detected
VlnPlot(PigE14.5C.8SOM, features = "nFeature_RNA", pt.size = 0.1) + NoLegend() #Bimodal distribution the majority of cells have around 800 or 4000 genes detected
VlnPlot(PigE14.5D.8SOM, features = "nFeature_RNA", pt.size = 0.1) + NoLegend() #Bimodal distribution the majority of cells have around 800 or 4000 genes detected
VlnPlot(PigE15A.10SOM, features = "nFeature_RNA", pt.size = 0.1) + NoLegend() #Bimodal distribution the majority of cells have around 800 or 4000 genes detected
VlnPlot(PigE15B.10SOM, features = "nFeature_RNA", pt.size = 0.1) + NoLegend() #Bimodal distribution the majority of cells have around 800 or 4000 genes detected
VlnPlot(PigE15C.10SOM, features = "nFeature_RNA", pt.size = 0.1) + NoLegend() #Bimodal distribution the majority of cells have around 800 or 4000 genes detected

#Identify and remove doublets with doublet finder
library(DoubletFinder)
#calculating pK values for each sample for doublet finder 
# PigE11.5A.EPS Pk 0.005
PigE11.5A.EPS_dbl <- paramSweep_v3(PigE11.5A.EPS, PCs = 1:50, sct = FALSE)
PigE11.5A.EPS_dbl_gt.calls <- summarizeSweep(PigE11.5A.EPS_dbl, GT = FALSE)
PigE11.5A.EPS.dblts.stats <- find.pK(PigE11.5A.EPS_dbl_gt.calls)
pK <- as.numeric(as.character(PigE11.5A.EPS.dblts.stats$pK))
BCmetric <- PigE11.5A.EPS.dblts.stats$BCmetric
pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
PigE11.5A.EPS_dbl_plot <- plot(x = pK, y = BCmetric, pch = 16,type="b",
                               col = "black",lty=1)+ abline(v=pK_choose,lwd=2,col='red',lty=2) +title("PigE11.5A.EPS") + text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.005*nrow(PigE11.5A.EPS@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
PigE11.5A.EPS <- doubletFinder_v3(PigE11.5A.EPS, PCs = 1:50, pN = 0.25, pK = 0.005, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(PigE11.5A.EPS, group.by= "DF.classifications_0.25_0.005_38")
#Set the Doublet finder as the default assay
PigE11.5A.EPS <- SetIdent(PigE11.5A.EPS, value = "DF.classifications_0.25_0.005_38")
#Keep singlets only
PigE11.5A.EPS <- subset(PigE11.5A.EPS, idents = c("Singlet"))
DimPlot(PigE11.5A.EPS, group.by= "DF.classifications_0.25_0.005_38")
FeaturePlot(PigE11.5A.EPS, features = "pANN_0.25_0.005_38")

# PigE11.5B.EPS Pk 0.005-------------------------------------------------------------------------------------
PigE11.5B.EPS_dbl <- paramSweep_v3(PigE11.5B.EPS, PCs = 1:50, sct = FALSE)
PigE11.5B.EPS_dbl_gt.calls <- summarizeSweep(PigE11.5B.EPS_dbl, GT = FALSE)
PigE11.5B.EPS.dblts.stats <- find.pK(PigE11.5B.EPS_dbl_gt.calls)
pK <- as.numeric(as.character(PigE11.5B.EPS.dblts.stats$pK))
BCmetric <- PigE11.5B.EPS.dblts.stats$BCmetric
pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
PigE11.5B.EPS_dbl_plot <- plot(x = pK, y = BCmetric, pch = 16,type="b",
                               col = "blue",lty=1)+ abline(v=pK_choose,lwd=2,col='red',lty=2) +title("The BCmvn distributions") + text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
PigE11.5B.EPS_dbl_plot
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.005*nrow(PigE11.5B.EPS@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
PigE11.5B.EPS <- doubletFinder_v3(PigE11.5B.EPS, PCs = 1:50, pN = 0.25, pK = 0.005, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(PigE11.5B.EPS, group.by= "DF.classifications_0.25_0.005_38")
#Set the Doublet finder as the default assay
PigE11.5B.EPS <- SetIdent(PigE11.5B.EPS, value = "DF.classifications_0.25_0.005_38")
#Keep singlets only
PigE11.5B.EPS <- subset(PigE11.5B.EPS, idents = c("Singlet"))
DimPlot(PigE11.5B.EPS, group.by= "DF.classifications_0.25_0.005_38")

# PigE12A.PS.EPS Pk 0.24-------------------------------------------------------------------------------------
PigE12A.PS_dbl <- paramSweep_v3(PigE12A.PS, PCs = 1:50, sct = FALSE)
PigE12A.PS_dbl_gt.calls <- summarizeSweep(PigE12A.PS_dbl, GT = FALSE)
PigE12A.PS.dblts.stats <- find.pK(PigE12A.PS_dbl_gt.calls)
pK <- as.numeric(as.character(PigE12A.PS.dblts.stats$pK))
BCmetric <- PigE12A.PS.dblts.stats$BCmetric
pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
PigE12A.PS_dbl_plot <- plot(x = pK, y = BCmetric, pch = 16,type="b",
                            col = "blue",lty=1)+ abline(v=pK_choose,lwd=2,col='red',lty=2) +title("The BCmvn distributions") + text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
PigE12A.PS_dbl_plot
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.24*nrow(PigE12A.PS@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
PigE12A.PS <- doubletFinder_v3(PigE12A.PS, PCs = 1:50, pN = 0.25, pK = 0.24, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(PigE12A.PS, group.by= "DF.classifications_0.25_0.24_952")
#Set the Doublet finder as the default assay
PigE12A.PS <- SetIdent(PigE12A.PS, value = "DF.classifications_0.25_0.24_952")
#Keep singlets only
PigE12A.PS <- subset(PigE12A.PS, idents = c("Singlet"))
DimPlot(PigE12A.PS, group.by= "DF.classifications_0.25_0.24_952")
FeaturePlot(PigE12A.PS, features = "pANN_0.25_0.24_952")
PigE12A.PS@meta.data[["DF.classifications_0.25_0.09_375"]]

# PigE12B.PS.EPS Pk 0.30-------------------------------------------------------------------------------------
PigE12B.PS_dbl <- paramSweep_v3(PigE12B.PS , PCs = 1:50, sct = FALSE)
PigE12B.PS_dbl_gt.calls <- summarizeSweep(PigE12B.PS_dbl, GT = FALSE)
PigE12B.PS.dblts.stats <- find.pK(PigE12B.PS_dbl_gt.calls)
pK <- as.numeric(as.character(PigE12B.PS.dblts.stats$pK))
BCmetric <- PigE12B.PS.dblts.stats$BCmetric
pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
PigE12B.PS_dbl_plot <- plot(x = pK, y = BCmetric, pch = 16,type="b",
                            col = "blue",lty=1)+ abline(v=pK_choose,lwd=2,col='red',lty=2) +title("The BCmvn distributions") + text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
PigE12B.PS_dbl_plot
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(PigE12B.PS@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
PigE12B.PS <- doubletFinder_v3(PigE12B.PS, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(PigE12B.PS, group.by= "DF.classifications_0.25_0.09_376")
#Set the Doublet finder as the default assay
PigE12B.PS <- SetIdent(PigE12B.PS, value = "DF.classifications_0.25_0.09_376")
#Keep singlets only
PigE12B.PS <- subset(PigE12B.PS, idents = c("Singlet"))
DimPlot(PigE12B.PS, group.by= "DF.classifications_0.25_0.09_376")

# PigE12.5A.PS2 Pk 0.08-------------------------------------------------------------------------------------
PigE12.5A.PS2_dbl <- paramSweep_v3(PigE12.5A.PS2, PCs = 1:50, sct = FALSE)
PigE12.5A.PS2_dbl_gt.calls <- summarizeSweep(PigE12.5A.PS2_dbl, GT = FALSE)
PigE12.5A.PS2.dblts.stats <- find.pK(PigE12.5A.PS2_dbl_gt.calls)
pK <- as.numeric(as.character(PigE12.5A.PS2.dblts.stats$pK))
BCmetric <- PigE12.5A.PS2.dblts.stats$BCmetric
pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
PigE12.5A.PS2_dbl_plot <- plot(x = pK, y = BCmetric, pch = 16,type="b",
                               col = "blue",lty=1)+ abline(v=pK_choose,lwd=2,col='red',lty=2) +title("The BCmvn distributions") + text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
PigE12.5A.PS2_dbl_plot
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.08*nrow(PigE12.5A.PS2@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
PigE12.5A.PS2 <- doubletFinder_v3(PigE12.5A.PS2, PCs = 1:50, pN = 0.25, pK = 0.08, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(PigE12.5A.PS2, group.by= "DF.classifications_0.25_0.08_147")
#Set the Doublet finder as the default assay
PigE12.5A.PS2 <- SetIdent(PigE12.5A.PS2, value = "DF.classifications_0.25_0.08_147")
#Keep singlets only
PigE12.5A.PS2 <- subset(PigE12.5A.PS2, idents = c("Singlet"))
DimPlot(PigE12.5A.PS2, group.by= "DF.classifications_0.25_0.08_147")

# PigE12.5B.PS2 Pk 0.05-------------------------------------------------------------------------------------
PigE12.5B.PS2_dbl <- paramSweep_v3(PigE12.5B.PS2, PCs = 1:50, sct = FALSE)
PigE12.5B.PS2_dbl_gt.calls <- summarizeSweep(PigE12.5B.PS2_dbl, GT = FALSE)
PigE12.5B.PS2.dblts.stats <- find.pK(PigE12.5B.PS2_dbl_gt.calls)
pK <- as.numeric(as.character(PigE12.5B.PS2.dblts.stats$pK))
BCmetric <- PigE12.5B.PS2.dblts.stats$BCmetric
pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
PigE12.5B.PS2_dbl_plot <- plot(x = pK, y = BCmetric, pch = 16,type="b",
                               col = "blue",lty=1)+ abline(v=pK_choose,lwd=2,col='red',lty=2) +title("The BCmvn distributions") + text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
PigE12.5B.PS2_dbl_plot
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.05*nrow(PigE12.5B.PS2@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
PigE12.5B.PS2 <- doubletFinder_v3(PigE12.5B.PS2, PCs = 1:50, pN = 0.25, pK = 0.05, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(PigE12.5B.PS2, group.by= "DF.classifications_0.25_0.05_61")
#Set the Doublet finder as the default assay
PigE12.5B.PS2 <- SetIdent(PigE12.5B.PS2, value = "DF.classifications_0.25_0.05_61")
#Keep singlets only
PigE12.5B.PS2 <- subset(PigE12.5B.PS2, idents = c("Singlet"))
DimPlot(PigE12.5B.PS2, group.by= "DF.classifications_0.25_0.05_61")

# PigE13A.NP Pk 0.07-------------------------------------------------------------------------------------
PigE13A.NP_dbl <- paramSweep_v3(PigE13A.NP, PCs = 1:50, sct = FALSE)
PigE13A.NP_dbl_gt.calls <- summarizeSweep(PigE13A.NP_dbl, GT = FALSE)
PigE13A.NP.dblts.stats <- find.pK(PigE13A.NP_dbl_gt.calls)
pK <- as.numeric(as.character(PigE13A.NP.dblts.stats$pK))
BCmetric <- PigE13A.NP.dblts.stats$BCmetric
pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
PigE13A.NP_dbl_plot <- plot(x = pK, y = BCmetric, pch = 16,type="b",
                            col = "blue",lty=1)+ abline(v=pK_choose,lwd=2,col='red',lty=2) +title("The BCmvn distributions") + text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
PigE13A.NP_dbl_plot 
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.07*nrow(PigE13A.NP@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
PigE13A.NP <- doubletFinder_v3(PigE13A.NP, PCs = 1:50, pN = 0.25, pK = 0.07, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(PigE13A.NP, group.by= "DF.classifications_0.25_0.07_57")
#Set the Doublet finder as the default assay
PigE13A.NP <- SetIdent(PigE13A.NP, value = "DF.classifications_0.25_0.07_57")
#Keep singlets only
PigE13A.NP <- subset(PigE13A.NP, idents = c("Singlet"))
DimPlot(PigE13A.NP, group.by= "DF.classifications_0.25_0.07_57")

# PigE13B.NP Pk 0.09-------------------------------------------------------------------------------------
PigE13B.NP_dbl <- paramSweep_v3(PigE13B.NP, PCs = 1:50, sct = FALSE)
PigE13B.NP_dbl_gt.calls <- summarizeSweep(PigE13B.NP_dbl, GT = FALSE)
PigE13B.NP.dblts.stats <- find.pK(PigE13B.NP_dbl_gt.calls)
pK <- as.numeric(as.character(PigE13B.NP.dblts.stats$pK))
BCmetric <- PigE13B.NP.dblts.stats$BCmetric
pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
PigE13B.NP_dbl_plot <- plot(x = pK, y = BCmetric, pch = 16,type="b",
                            col = "blue",lty=1)+ abline(v=pK_choose,lwd=2,col='red',lty=2) +title("The BCmvn distributions") + text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
PigE13B.NP_dbl_plot
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.09*nrow(PigE13B.NP@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
PigE13B.NP <- doubletFinder_v3(PigE13B.NP, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(PigE13B.NP, group.by= "DF.classifications_0.25_0.09_557")
#Set the Doublet finder as the default assay
PigE13B.NP <- SetIdent(PigE13B.NP, value = "DF.classifications_0.25_0.09_557")
#Keep singlets only
PigE13B.NP <- subset(PigE13B.NP, idents = c("Singlet"))
DimPlot(PigE13B.NP, group.by= "DF.classifications_0.25_0.09_557")

# PigE13C.NP Pk 0.1-------------------------------------------------------------------------------------
PigE13C.NP_dbl <- paramSweep_v3(PigE13C.NP, PCs = 1:50, sct = FALSE)
PigE13C.NP_dbl_gt.calls <- summarizeSweep(PigE13C.NP_dbl, GT = FALSE)
PigE13C.NP.dblts.stats <- find.pK(PigE13C.NP_dbl_gt.calls)
pK <- as.numeric(as.character(PigE13C.NP.dblts.stats$pK))
BCmetric <- PigE13C.NP.dblts.stats$BCmetric
pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
PigE13C.NP_dbl_plot <- plot(x = pK, y = BCmetric, pch = 16,type="b",
                            col = "blue",lty=1)+ abline(v=pK_choose,lwd=2,col='red',lty=2) +title("The BCmvn distributions") + text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
PigE13C.NP_dbl_plot
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.1*nrow(PigE13C.NP@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
PigE13C.NP <- doubletFinder_v3(PigE13C.NP, PCs = 1:50, pN = 0.25, pK = 0.1, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(PigE13C.NP, group.by= "DF.classifications_0.25_0.1_550")
#Set the Doublet finder as the default assay
PigE13C.NP <- SetIdent(PigE13C.NP, value = "DF.classifications_0.25_0.1_550")
#Keep singlets only
PigE13C.NP <- subset(PigE13C.NP, idents = c("Singlet"))
DimPlot(PigE13C.NP, group.by= "DF.classifications_0.25_0.1_550")

# PigE13D.NP Pk 0.16-------------------------------------------------------------------------------------
PigE13D.NP_dbl <- paramSweep_v3(PigE13D.NP, PCs = 1:50, sct = FALSE)
PigE13D.NP_dbl_gt.calls <- summarizeSweep(PigE13D.NP_dbl, GT = FALSE)
PigE13D.NP.dblts.stats <- find.pK(PigE13D.NP_dbl_gt.calls)
pK <- as.numeric(as.character(PigE13D.NP.dblts.stats$pK))
BCmetric <- PigE13D.NP.dblts.stats$BCmetric
pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
PigE13D.NP_dbl_plot <- plot(x = pK, y = BCmetric, pch = 16,type="b",
                            col = "blue",lty=1)+ abline(v=pK_choose,lwd=2,col='red',lty=2) +title("The BCmvn distributions") + text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
PigE13D.NP_dbl_plot
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.16*nrow(PigE13D.NP@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
PigE13D.NP <- doubletFinder_v3(PigE13D.NP, PCs = 1:50, pN = 0.25, pK = 0.16, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(PigE13D.NP, group.by= "DF.classifications_0.25_0.16_3641")
#Set the Doublet finder as the default assay
PigE13D.NP <- SetIdent(PigE13D.NP, value = "DF.classifications_0.25_0.16_3641")
#Keep singlets only
PigE13D.NP <- subset(PigE13D.NP, idents = c("Singlet"))
DimPlot(PigE13D.NP, group.by= "DF.classifications_0.25_0.16_3641")
FeaturePlot(PigE13D.NP, features = "pANN_0.25_0.09_375")
PigE13D.NP@meta.data[["DF.classifications_0.25_0.09_375"]]

# PigE13E.NP Pk 0.15-------------------------------------------------------------------------------------
PigE13E.NP_dbl <- paramSweep_v3(PigE13E.NP, PCs = 1:50, sct = FALSE)
PigE13E.NP_dbl_gt.calls <- summarizeSweep(PigE13E.NP_dbl, GT = FALSE)
PigE13E.NP.dblts.stats <- find.pK(PigE13E.NP_dbl_gt.calls)
pK <- as.numeric(as.character(PigE13E.NP.dblts.stats$pK))
BCmetric <- PigE13E.NP.dblts.stats$BCmetric
pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
PigE13E.NP_dbl_plot <- plot(x = pK, y = BCmetric, pch = 16,type="b",
                            col = "blue",lty=1)+ abline(v=pK_choose,lwd=2,col='red',lty=2) +title("The BCmvn distributions") + text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
PigE13E.NP_dbl_plot #UP TO HERE!!!!!!!!!!!
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.15*nrow(PigE13E.NP@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
PigE13E.NP <- doubletFinder_v3(PigE13E.NP, PCs = 1:50, pN = 0.25, pK = 0.15, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(PigE13E.NP, group.by= "DF.classifications_0.25_0.15_3125")
#Set the Doublet finder as the default assay
PigE13E.NP <- SetIdent(PigE13E.NP, value = "DF.classifications_0.25_0.15_3125")
#Keep singlets only
PigE13E.NP <- subset(PigE13E.NP, idents = c("Singlet"))
DimPlot(PigE13E.NP, group.by= "DF.classifications_0.25_0.15_3125")

# PigE13.5A.NG Pk 0.2-------------------------------------------------------------------------------------
PigE13.5A.NG_dbl <- paramSweep_v3(PigE13.5A.NG, PCs = 1:50, sct = FALSE)
PigE13.5A.NG_dbl_gt.calls <- summarizeSweep(PigE13.5A.NG_dbl, GT = FALSE)
PigE13.5A.NG.dblts.stats <- find.pK(PigE13.5A.NG_dbl_gt.calls)
pK <- as.numeric(as.character(PigE13.5A.NG.dblts.stats$pK))
BCmetric <- PigE13.5A.NG.dblts.stats$BCmetric
pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
PigE13.5A.NG_dbl_plot <- plot(x = pK, y = BCmetric, pch = 16,type="b",
                              col = "blue",lty=1)+ abline(v=pK_choose,lwd=2,col='red',lty=2) +title("The BCmvn distributions") + text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
PigE13.5A.NG_dbl_plot
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.2*nrow(PigE13.5A.NG@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
PigE13.5A.NG <- doubletFinder_v3(PigE13.5A.NG, PCs = 1:50, pN = 0.25, pK = 0.2, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(PigE13.5A.NG, group.by= "DF.classifications_0.25_0.2_768")
#Set the Doublet finder as the default assay
PigE13.5A.NG <- SetIdent(PigE13.5A.NG, value = "DF.classifications_0.25_0.2_768")
#Keep singlets only
PigE13.5A.NG <- subset(PigE13.5A.NG, idents = c("Singlet"))
DimPlot(PigE13.5A.NG, group.by= "DF.classifications_0.25_0.2_768")

# PigE13.5B.NG Pk 0.24-------------------------------------------------------------------------------------
PigE13.5B.NG_dbl <- paramSweep_v3(PigE13.5B.NG, PCs = 1:50, sct = FALSE)
PigE13.5B.NG_dbl_gt.calls <- summarizeSweep(PigE13.5B.NG_dbl, GT = FALSE)
PigE13.5B.NG.dblts.stats <- find.pK(PigE13.5B.NG_dbl_gt.calls)
pK <- as.numeric(as.character(PigE13.5B.NG.dblts.stats$pK))
BCmetric <- PigE13.5B.NG.dblts.stats$BCmetric
pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
PigE13.5B.NG_dbl_plot <- plot(x = pK, y = BCmetric, pch = 16,type="b",
                              col = "blue",lty=1)+ abline(v=pK_choose,lwd=2,col='red',lty=2) +title("The BCmvn distributions") + text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
PigE13.5B.NG_dbl_plot
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.24*nrow(PigE13.5B.NG@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
PigE13.5B.NG <- doubletFinder_v3(PigE13.5B.NG, PCs = 1:50, pN = 0.25, pK = 0.24, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(PigE13.5B.NG, group.by= "DF.classifications_0.25_0.24_826")
#Set the Doublet finder as the default assay
PigE13.5B.NG <- SetIdent(PigE13.5B.NG, value = "DF.classifications_0.25_0.24_826")
#Keep singlets only
PigE13.5B.NG <- subset(PigE13.5B.NG, idents = c("Singlet"))
DimPlot(PigE13.5B.NG, group.by= "DF.classifications_0.25_0.24_826")

# PigE13.5C.NG Pk 0.14-------------------------------------------------------------------------------------
PigE13.5C.NG_dbl <- paramSweep_v3(PigE13.5C.NG, PCs = 1:50, sct = FALSE)
PigE13.5C.NG_dbl_gt.calls <- summarizeSweep(PigE13.5C.NG_dbl, GT = FALSE)
PigE13.5C.NG.dblts.stats <- find.pK(PigE13.5C.NG_dbl_gt.calls)
pK <- as.numeric(as.character(PigE13.5C.NG.dblts.stats$pK))
BCmetric <- PigE13.5C.NG.dblts.stats$BCmetric
pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
PigE13.5C.NG_dbl_plot <- plot(x = pK, y = BCmetric, pch = 16,type="b",
                              col = "blue",lty=1)+ abline(v=pK_choose,lwd=2,col='red',lty=2) +title("The BCmvn distributions") + text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
PigE13.5C.NG_dbl_plot
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.14*nrow(PigE13.5C.NG@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
PigE13.5C.NG <- doubletFinder_v3(PigE13.5C.NG, PCs = 1:50, pN = 0.25, pK = 0.14, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(PigE13.5C.NG, group.by= "DF.classifications_0.25_0.14_1188")
#Set the Doublet finder as the default assay
PigE13.5C.NG <- SetIdent(PigE13.5C.NG, value = "DF.classifications_0.25_0.14_1188")
#Keep singlets only
PigE13.5C.NG <- subset(PigE13.5C.NG, idents = c("Singlet"))
DimPlot(PigE13.5C.NG, group.by= "DF.classifications_0.25_0.14_1188")

# PigE13.5D.NG Pk 0.23 -------------------------------------------------------------------------------------
PigE13.5D.NG_dbl <- paramSweep_v3(PigE13.5D.NG, PCs = 1:50, sct = FALSE)
PigE13.5D.NG_dbl_gt.calls <- summarizeSweep(PigE13.5D.NG_dbl, GT = FALSE)
PigE13.5D.NG.dblts.stats <- find.pK(PigE13.5D.NG_dbl_gt.calls)
pK <- as.numeric(as.character(PigE13.5D.NG.dblts.stats$pK))
BCmetric <- PigE13.5D.NG.dblts.stats$BCmetric
pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
PigE13.5D.NG_dbl_plot <- plot(x = pK, y = BCmetric, pch = 16,type="b",
                              col = "blue",lty=1)+ abline(v=pK_choose,lwd=2,col='red',lty=2) +title("The BCmvn distributions") + text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
PigE13.5D.NG_dbl_plot
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.23*nrow(PigE13.5D.NG@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
PigE13.5D.NG <- doubletFinder_v3(PigE13.5D.NG, PCs = 1:50, pN = 0.25, pK = 0.23, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(PigE13.5D.NG, group.by= "DF.classifications_0.25_0.23_2805")
#Set the Doublet finder as the default assay
PigE13.5D.NG <- SetIdent(PigE13.5D.NG, value = "DF.classifications_0.25_0.23_2805")
#Keep singlets only
PigE13.5D.NG <- subset(PigE13.5D.NG, idents = c("Singlet"))
DimPlot(PigE13.5D.NG, group.by= "DF.classifications_0.25_0.23_2805")

# PigE14A.6SOM Pk 0.28-------------------------------------------------------------------------------------
PigE14A.6SOM_dbl <- paramSweep_v3(PigE14A.6SOM, PCs = 1:50, sct = FALSE)
PigE14A.6SOM_dbl_gt.calls <- summarizeSweep(PigE14A.6SOM_dbl, GT = FALSE)
PigE14A.6SOM.dblts.stats <- find.pK(PigE14A.6SOM_dbl_gt.calls)
pK <- as.numeric(as.character(PigE14A.6SOM.dblts.stats$pK))
BCmetric <- PigE14A.6SOM.dblts.stats$BCmetric
pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
PigE14A.6SOM_dbl_plot <- plot(x = pK, y = BCmetric, pch = 16,type="b",
                              col = "blue",lty=1)+ abline(v=pK_choose,lwd=2,col='red',lty=2) +title("The BCmvn distributions") + text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
PigE14A.6SOM_dbl_plot
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.28*nrow(PigE14A.6SOM@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
PigE14A.6SOM <- doubletFinder_v3(PigE14A.6SOM, PCs = 1:50, pN = 0.25, pK = 0.28, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(PigE14A.6SOM, group.by= "DF.classifications_0.25_0.28_661")
#Set the Doublet finder as the default assay
PigE14A.6SOM <- SetIdent(PigE14A.6SOM, value = "DF.classifications_0.25_0.28_661")
#Keep singlets only
PigE14A.6SOM <- subset(PigE14A.6SOM, idents = c("Singlet"))
DimPlot(PigE14A.6SOM, group.by= "DF.classifications_0.25_0.28_661")

# PigE14B.6SOM Pk 0.28-------------------------------------------------------------------------------------
PigE14B.6SOM_dbl <- paramSweep_v3(PigE14B.6SOM, PCs = 1:50, sct = FALSE)
PigE14B.6SOM_dbl_gt.calls <- summarizeSweep(PigE14B.6SOM_dbl, GT = FALSE)
PigE14B.6SOM.dblts.stats <- find.pK(PigE14B.6SOM_dbl_gt.calls)
pK <- as.numeric(as.character(PigE14B.6SOM.dblts.stats$pK))
BCmetric <- PigE14B.6SOM.dblts.stats$BCmetric
pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
PigE14B.6SOM_dbl_plot <- plot(x = pK, y = BCmetric, pch = 16,type="b",
                              col = "blue",lty=1)+ abline(v=pK_choose,lwd=2,col='red',lty=2) +title("The BCmvn distributions") + text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
PigE14B.6SOM_dbl_plot
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.28*nrow(PigE14B.6SOM@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
PigE14B.6SOM <- doubletFinder_v3(PigE14B.6SOM, PCs = 1:50, pN = 0.25, pK = 0.28, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(PigE14B.6SOM, group.by= "DF.classifications_0.25_0.28_588")
#Set the Doublet finder as the default assay
PigE14B.6SOM <- SetIdent(PigE14B.6SOM, value = "DF.classifications_0.25_0.28_588")
#Keep singlets only
PigE14B.6SOM <- subset(PigE14B.6SOM, idents = c("Singlet"))
DimPlot(PigE14B.6SOM, group.by= "DF.classifications_0.25_0.28_588")

# PigE14.5A.8SOM Pk 0.03-------------------------------------------------------------------------------------
PigE14.5A.8SOM_dbl <- paramSweep_v3(PigE14.5A.8SOM, PCs = 1:50, sct = FALSE)
PigE14.5A.8SOM_dbl_gt.calls <- summarizeSweep(PigE14.5A.8SOM_dbl, GT = FALSE)
PigE14.5A.8SOM.dblts.stats <- find.pK(PigE14.5A.8SOM_dbl_gt.calls)
pK <- as.numeric(as.character(PigE14.5A.8SOM.dblts.stats$pK))
BCmetric <- PigE14.5A.8SOM.dblts.stats$BCmetric
pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
PigE14.5A.8SOM_dbl_plot <- plot(x = pK, y = BCmetric, pch = 16,type="b",
                                col = "blue",lty=1)+ abline(v=pK_choose,lwd=2,col='red',lty=2) +title("The BCmvn distributions") + text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
PigE14.5A.8SOM_dbl_plot
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.03*nrow(PigE14.5A.8SOM@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
PigE14.5A.8SOM <- doubletFinder_v3(PigE14.5A.8SOM, PCs = 1:50, pN = 0.25, pK = 0.03, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(PigE14.5A.8SOM, group.by= "DF.classifications_0.25_0.03_120")
#Set the Doublet finder as the default assay
PigE14.5A.8SOM <- SetIdent(PigE14.5A.8SOM, value = "DF.classifications_0.25_0.03_120")
#Keep singlets only
PigE14.5A.8SOM <- subset(PigE14.5A.8SOM, idents = c("Singlet"))
DimPlot(PigE14.5A.8SOM, group.by= "DF.classifications_0.25_0.03_120")

# PigE14.5B.8SOM Pk 0.03-------------------------------------------------------------------------------------
PigE14.5B.8SOM_dbl <- paramSweep_v3(PigE14.5B.8SOM, PCs = 1:50, sct = FALSE)
PigE14.5B.8SOM_dbl_gt.calls <- summarizeSweep(PigE14.5B.8SOM_dbl, GT = FALSE)
PigE14.5B.8SOM.dblts.stats <- find.pK(PigE14.5B.8SOM_dbl_gt.calls)
pK <- as.numeric(as.character(PigE14.5B.8SOM.dblts.stats$pK))
BCmetric <- PigE14.5B.8SOM.dblts.stats$BCmetric
pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
PigE14.5B.8SOM_dbl_plot <- plot(x = pK, y = BCmetric, pch = 16,type="b",
                                col = "blue",lty=1)+ abline(v=pK_choose,lwd=2,col='red',lty=2) +title("The BCmvn distributions") + text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
PigE14.5B.8SOM_dbl_plot
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.03*nrow(PigE14.5B.8SOM@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
PigE14.5B.8SOM <- doubletFinder_v3(PigE14.5B.8SOM, PCs = 1:50, pN = 0.25, pK = 0.03, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(PigE14.5B.8SOM, group.by= "DF.classifications_0.25_0.03_61")
#Set the Doublet finder as the default assay
PigE14.5B.8SOM <- SetIdent(PigE14.5B.8SOM, value = "DF.classifications_0.25_0.03_61")
#Keep singlets only
PigE14.5B.8SOM <- subset(PigE14.5B.8SOM, idents = c("Singlet"))
DimPlot(PigE14.5B.8SOM, group.by= "DF.classifications_0.25_0.03_61")

# PigE14.5C.8SOM Pk 0.27-------------------------------------------------------------------------------------
PigE14.5C.8SOM_dbl <- paramSweep_v3(PigE14.5C.8SOM, PCs = 1:50, sct = FALSE)
PigE14.5C.8SOM_dbl_gt.calls <- summarizeSweep(PigE14.5C.8SOM_dbl, GT = FALSE)
PigE14.5C.8SOM.dblts.stats <- find.pK(PigE14.5C.8SOM_dbl_gt.calls)
pK <- as.numeric(as.character(PigE14.5C.8SOM.dblts.stats$pK))
BCmetric <- PigE14.5C.8SOM.dblts.stats$BCmetric
pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
PigE14.5C.8SOM_dbl_plot <- plot(x = pK, y = BCmetric, pch = 16,type="b",
                                col = "blue",lty=1)+ abline(v=pK_choose,lwd=2,col='red',lty=2) +title("The BCmvn distributions") + text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
PigE14.5C.8SOM_dbl_plot
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.27*nrow(PigE14.5C.8SOM@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
PigE14.5C.8SOM <- doubletFinder_v3(PigE14.5C.8SOM, PCs = 1:50, pN = 0.25, pK = 0.27, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(PigE14.5C.8SOM, group.by= "DF.classifications_0.25_0.27_374")
#Set the Doublet finder as the default assay
PigE14.5C.8SOM <- SetIdent(PigE14.5C.8SOM, value = "DF.classifications_0.25_0.27_374")
#Keep singlets only
PigE14.5C.8SOM <- subset(PigE14.5C.8SOM, idents = c("Singlet"))
DimPlot(PigE14.5C.8SOM, group.by= "DF.classifications_0.25_0.27_374")

# PigE14.5D.8SOM Pk 0.30-------------------------------------------------------------------------------------
PigE14.5D.8SOM_dbl <- paramSweep_v3(PigE14.5D.8SOM, PCs = 1:50, sct = FALSE)
PigE14.5D.8SOM_dbl_gt.calls <- summarizeSweep(PigE14.5D.8SOM_dbl, GT = FALSE)
PigE14.5D.8SOM.dblts.stats <- find.pK(PigE14.5D.8SOM_dbl_gt.calls)
pK <- as.numeric(as.character(PigE14.5D.8SOM.dblts.stats$pK))
BCmetric <- PigE14.5D.8SOM.dblts.stats$BCmetric
pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
PigE14.5D.8SOM_dbl_plot <- plot(x = pK, y = BCmetric, pch = 16,type="b",
                                col = "blue",lty=1)+ abline(v=pK_choose,lwd=2,col='red',lty=2) +title("The BCmvn distributions") + text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
PigE14.5D.8SOM_dbl_plot
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.30*nrow(PigE14.5D.8SOM@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
PigE14.5D.8SOM <- doubletFinder_v3(PigE14.5D.8SOM, PCs = 1:50, pN = 0.25, pK = 0.30, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(PigE14.5D.8SOM, group.by= "DF.classifications_0.25_0.3_420")
#Set the Doublet finder as the default assay
PigE14.5D.8SOM <- SetIdent(PigE14.5D.8SOM, value = "DF.classifications_0.25_0.3_420")
#Keep singlets only
PigE14.5D.8SOM <- subset(PigE14.5D.8SOM, idents = c("Singlet"))
DimPlot(PigE14.5D.8SOM, group.by= "DF.classifications_0.25_0.3_420")

# PigE15A.10SOM Pk 0.05-------------------------------------------------------------------------------------
PigE15A.10SOM_dbl <- paramSweep_v3(PigE15A.10SOM , PCs = 1:50, sct = FALSE)
PigE15A.10SOM_dbl_gt.calls <- summarizeSweep(PigE15A.10SOM_dbl, GT = FALSE)
PigE15A.10SOM.dblts.stats <- find.pK(PigE15A.10SOM_dbl_gt.calls)
pK <- as.numeric(as.character(PigE15A.10SOM.dblts.stats$pK))
BCmetric <- PigE15A.10SOM.dblts.stats$BCmetric
pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
PigE15A.10SOM_dbl_plot <- plot(x = pK, y = BCmetric, pch = 16,type="b",
                               col = "blue",lty=1)+ abline(v=pK_choose,lwd=2,col='red',lty=2) +title("The BCmvn distributions") + text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
PigE15A.10SOM_dbl_plot
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.05*nrow(PigE15A.10SOM@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
PigE15A.10SOM <- doubletFinder_v3(PigE15A.10SOM, PCs = 1:50, pN = 0.25, pK = 0.05, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(PigE15A.10SOM, group.by= "DF.classifications_0.25_0.05_163")
#Set the Doublet finder as the default assay
PigE15A.10SOM <- SetIdent(PigE15A.10SOM, value = "DF.classifications_0.25_0.05_163")
#Keep singlets only
PigE15A.10SOM <- subset(PigE15A.10SOM, idents = c("Singlet"))
DimPlot(PigE15A.10SOM, group.by= "DF.classifications_0.25_0.05_163")

# PigE15B.10SOM Pk 0.04-------------------------------------------------------------------------------------
PigE15B.10SOM_dbl <- paramSweep_v3(PigE15B.10SOM, PCs = 1:50, sct = FALSE)
PigE15B.10SOM_dbl_gt.calls <- summarizeSweep(PigE15B.10SOM_dbl, GT = FALSE)
PigE15B.10SOM.dblts.stats <- find.pK(PigE15B.10SOM_dbl_gt.calls)
pK <- as.numeric(as.character(PigE15B.10SOM.dblts.stats$pK))
BCmetric <- PigE15B.10SOM.dblts.stats$BCmetric
pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
PigE15B.10SOM_dbl_plot <- plot(x = pK, y = BCmetric, pch = 16,type="b",
                               col = "blue",lty=1)+ abline(v=pK_choose,lwd=2,col='red',lty=2) +title("The BCmvn distributions") + text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
PigE15B.10SOM_dbl_plot
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.04*nrow(PigE15B.10SOM@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
PigE15B.10SOM <- doubletFinder_v3(PigE15B.10SOM, PCs = 1:50, pN = 0.25, pK = 0.04, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(PigE15B.10SOM, group.by= "DF.classifications_0.25_0.04_211")
#Set the Doublet finder as the default assay
PigE15B.10SOM <- SetIdent(PigE15B.10SOM, value = "DF.classifications_0.25_0.04_211")
#Keep singlets only
PigE15B.10SOM <- subset(PigE15B.10SOM, idents = c("Singlet"))
DimPlot(PigE15B.10SOM, group.by= "DF.classifications_0.25_0.04_211")

# PigE15C.10SOM Pk 0.24-------------------------------------------------------------------------------------
PigE15C.10SOM_dbl <- paramSweep_v3(PigE15C.10SOM, PCs = 1:50, sct = FALSE)
PigE15C.10SOM_dbl_gt.calls <- summarizeSweep(PigE15C.10SOM_dbl, GT = FALSE)
PigE15C.10SOM.dblts.stats <- find.pK(PigE15C.10SOM_dbl_gt.calls)
pK <- as.numeric(as.character(PigE15C.10SOM.dblts.stats$pK))
BCmetric <- PigE15C.10SOM.dblts.stats$BCmetric
pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
PigE15C.10SOM_dbl_plot <- plot(x = pK, y = BCmetric, pch = 16,type="b",
                               col = "blue",lty=1)+ abline(v=pK_choose,lwd=2,col='red',lty=2) +title("The BCmvn distributions") + text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
PigE15C.10SOM_dbl_plot
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.24*nrow(PigE15C.10SOM@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
PigE15C.10SOM <- doubletFinder_v3(PigE15C.10SOM, PCs = 1:50, pN = 0.25, pK = 0.24, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(PigE15C.10SOM, group.by= "DF.classifications_0.25_0.24_1134")
#Set the Doublet finder as the default assay
PigE15C.10SOM <- SetIdent(PigE15C.10SOM, value = "DF.classifications_0.25_0.24_1134")
#Keep singlets only
PigE15C.10SOM <- subset(PigE15C.10SOM, idents = c("Singlet"))
DimPlot(PigE15C.10SOM, group.by= "DF.classifications_0.25_0.24_1134")

#Add a new column in the metadata to keep the original sample info (NOTE this is also stored in orig.ident but this column would ot be present in non-seurat datasets)
PigE11.5A.EPS$Sample <- "E11.5A.EPS"
PigE11.5B.EPS$Sample <- "E11.5B.EPS"
PigE12A.PS$Sample <-"E12A.PS" 
PigE12B.PS$Sample <-"E12B.PS" 
PigE12.5A.PS2$Sample <-"E12.5A.PS2" 
PigE12.5B.PS2$Sample <-"E12.5B.PS2" 
PigE13A.NP$Sample <-"E13A.NP" 
PigE13B.NP$Sample <-"E13B.NP" 
PigE13C.NP$Sample <-"E13C.NP" 
PigE13D.NP$Sample <-"E13D.NP" 
PigE13E.NP$Sample <-"E13E.NP" 
PigE13.5A.NG$Sample <-"E13.5A.NG" 
PigE13.5B.NG$Sample <-"E13.5B.NG" 
PigE13.5C.NG$Sample <-"E13.5C.NG"
PigE13.5D.NG$Sample <-"E13.5D.NG" 
PigE14A.6SOM$Sample <-"E14A.6SOM" 
PigE14B.6SOM$Sample <-"E14B.6SOM" 
PigE14.5A.8SOM$Sample <-"E14.5A.8SOM" 
PigE14.5B.8SOM$Sample <-"E14.5B.8SOM" 
PigE14.5C.8SOM$Sample <-"E14.5C.8SOM" 
PigE14.5D.8SOM$Sample <-"E14.5D.8SOM" 
PigE15A.10SOM$Sample <-"E15A.10SOM" 
PigE15B.10SOM$Sample <-"E15B.10SOM" 
PigE15C.10SOM$Sample <-"E15C.10SOM" 

# Add stage information to a column in the metadata
PigE11.5A.EPS$Stage <- "E11.5"
PigE11.5B.EPS$Stage <- "E11.5"
PigE12A.PS$Stage <-"E12" 
PigE12B.PS$Stage <-"E12" 
PigE12.5A.PS2$Stage <-"E12.5" 
PigE12.5B.PS2$Stage <-"E12.5" 
PigE13A.NP$Stage <-"E13" 
PigE13B.NP$Stage <-"E13" 
PigE13C.NP$Stage <-"E13" 
PigE13D.NP$Stage <-"E13" 
PigE13E.NP$Stage <-"E13" 
PigE13.5A.NG$Stage <-"E13.5" 
PigE13.5B.NG$Stage <-"E13.5" 
PigE13.5C.NG$Stage <-"E13.5"
PigE13.5D.NG$Stage <-"E13.5" 
PigE14A.6SOM$Stage <-"E14" 
PigE14B.6SOM$Stage <-"E14" 
PigE14.5A.8SOM$Stage <-"E14.5" 
PigE14.5B.8SOM$Stage <-"E14.5" 
PigE14.5C.8SOM$Stage <-"E14.5" 
PigE14.5D.8SOM$Stage <-"E14.5" 
PigE15A.10SOM$Stage <-"E15" 
PigE15B.10SOM$Stage <-"E15" 
PigE15C.10SOM$Stage <-"E15" 


# This section contains repetitive code
# The code below may be a better choice for pre-processing

#object_list <- c("PigE11.5A.EPS", "PigE11.5B.EPS", "PigE12A.PS", "PigE12B.PS","PigE12.5A.PS2", "PigE12.5B.PS2", "PigE13A.NP", "PigE13B.NP", "PigE13C.NP", "PigE13D.NP", "PigE13E.NP", "PigE13.5A.NG","PigE13.5B.NG", "PigE13.5C.NG", "PigE13.5D.NG", "PigE14A.6SOM","PigE14B.6SOM", "PigE14.5A.8SOM", "PigE14.5B.8SOM", "PigE14.5C.8SOM","PigE14.5D.8SOM", "PigE15A.10SOM", "PigE15B.10SOM", "PigE15C.10SOM")
#object_list <- lapply(X = object_list, FUN = function(x) {
#x <- NormalizeData(x, verbose = FALSE)
#x <- FindVariableFeatures(x, verbose = FALSE)
#x <- ScaleData(x, features = rownames(x), verbose = FALSE)
#x <- RunPCA(x, npcs = 100, verbose = FALSE)
#})

## To normalise the data ##
PigE11.5A.EPS <- NormalizeData(PigE11.5A.EPS, normalization.method = "LogNormalize", scale.factor = 10000)
PigE11.5B.EPS <- NormalizeData(PigE11.5B.EPS, normalization.method = "LogNormalize", scale.factor = 10000)
PigE12A.PS <- NormalizeData(PigE12A.PS, normalization.method = "LogNormalize", scale.factor = 10000)
PigE12B.PS <- NormalizeData(PigE12B.PS, normalization.method = "LogNormalize", scale.factor = 10000)
PigE12.5A.PS2 <- NormalizeData(PigE12.5A.PS2, normalization.method = "LogNormalize", scale.factor = 10000)
PigE12.5B.PS2 <- NormalizeData(PigE12.5B.PS2, normalization.method = "LogNormalize", scale.factor = 10000)
PigE13A.NP <- NormalizeData(PigE13A.NP, normalization.method = "LogNormalize", scale.factor = 10000)
PigE13B.NP <- NormalizeData(PigE13B.NP, normalization.method = "LogNormalize", scale.factor = 10000)
PigE13C.NP <- NormalizeData(PigE13C.NP, normalization.method = "LogNormalize", scale.factor = 10000)
PigE13D.NP <- NormalizeData(PigE13D.NP, normalization.method = "LogNormalize", scale.factor = 10000)
PigE13E.NP <- NormalizeData(PigE13E.NP, normalization.method = "LogNormalize", scale.factor = 10000)
PigE13.5A.NG <- NormalizeData(PigE13.5A.NG, normalization.method = "LogNormalize", scale.factor = 10000)
PigE13.5B.NG <- NormalizeData(PigE13.5B.NG, normalization.method = "LogNormalize", scale.factor = 10000)
PigE13.5C.NG <- NormalizeData(PigE13.5C.NG, normalization.method = "LogNormalize", scale.factor = 10000)
PigE13.5D.NG <- NormalizeData(PigE13.5D.NG, normalization.method = "LogNormalize", scale.factor = 10000)
PigE14A.6SOM <- NormalizeData(PigE14A.6SOM, normalization.method = "LogNormalize", scale.factor = 10000)
PigE14B.6SOM <- NormalizeData(PigE14B.6SOM, normalization.method = "LogNormalize", scale.factor = 10000)
PigE14.5A.8SOM <- NormalizeData(PigE14.5A.8SOM, normalization.method = "LogNormalize", scale.factor = 10000)
PigE14.5B.8SOM <- NormalizeData(PigE14.5B.8SOM, normalization.method = "LogNormalize", scale.factor = 10000)
PigE14.5C.8SOM <- NormalizeData(PigE14.5C.8SOM, normalization.method = "LogNormalize", scale.factor = 10000)
PigE14.5D.8SOM <- NormalizeData(PigE14.5D.8SOM, normalization.method = "LogNormalize", scale.factor = 10000)
PigE15A.10SOM <- NormalizeData(PigE15A.10SOM, normalization.method = "LogNormalize", scale.factor = 10000)
PigE15B.10SOM <- NormalizeData(PigE15B.10SOM, normalization.method = "LogNormalize", scale.factor = 10000)
PigE15C.10SOM <- NormalizeData(PigE15C.10SOM, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
PigE11.5A.EPS <- FindVariableFeatures(object = PigE11.5A.EPS, selection.method = "vst", nfeatures = 5000) #was 4000
PigE11.5B.EPS <- FindVariableFeatures(object = PigE11.5B.EPS, selection.method = "vst", nfeatures = 5000)
PigE12A.PS <- FindVariableFeatures(object = PigE12A.PS, selection.method = "vst", nfeatures = 5000)
PigE12B.PS <- FindVariableFeatures(object = PigE12B.PS, selection.method = "vst", nfeatures = 5000)
PigE12.5A.PS2 <- FindVariableFeatures(object = PigE12.5A.PS2, selection.method = "vst", nfeatures = 5000)
PigE12.5B.PS2 <- FindVariableFeatures(object = PigE12.5B.PS2, selection.method = "vst", nfeatures = 5000)
PigE13A.NP <- FindVariableFeatures(object = PigE13A.NP, selection.method = "vst", nfeatures = 5000)
PigE13B.NP <- FindVariableFeatures(object = PigE13B.NP, selection.method = "vst", nfeatures = 5000)
PigE13C.NP <- FindVariableFeatures(object = PigE13C.NP, selection.method = "vst", nfeatures = 5000)
PigE13D.NP <- FindVariableFeatures(object = PigE13D.NP, selection.method = "vst", nfeatures = 5000)
PigE13E.NP <- FindVariableFeatures(object = PigE13E.NP, selection.method = "vst", nfeatures = 5000)
PigE13.5A.NG <- FindVariableFeatures(object = PigE13.5A.NG, selection.method = "vst", nfeatures = 5000)
PigE13.5B.NG <- FindVariableFeatures(object = PigE13.5B.NG, selection.method = "vst", nfeatures = 5000)
PigE13.5C.NG <- FindVariableFeatures(object = PigE13.5C.NG, selection.method = "vst", nfeatures = 5000)
PigE13.5D.NG <- FindVariableFeatures(object = PigE13.5D.NG, selection.method = "vst", nfeatures = 5000)
PigE14A.6SOM <- FindVariableFeatures(object = PigE14A.6SOM, selection.method = "vst", nfeatures = 5000)
PigE14B.6SOM <- FindVariableFeatures(object = PigE14B.6SOM, selection.method = "vst", nfeatures = 5000)
PigE14.5A.8SOM <- FindVariableFeatures(object = PigE14.5A.8SOM, selection.method = "vst", nfeatures = 5000)
PigE14.5B.8SOM <- FindVariableFeatures(object = PigE14.5B.8SOM, selection.method = "vst", nfeatures = 5000)
PigE14.5C.8SOM <- FindVariableFeatures(object = PigE14.5C.8SOM, selection.method = "vst", nfeatures = 5000)
PigE14.5D.8SOM <- FindVariableFeatures(object = PigE14.5D.8SOM, selection.method = "vst", nfeatures = 5000)
PigE15A.10SOM <- FindVariableFeatures(object = PigE15A.10SOM, selection.method = "vst", nfeatures = 5000)
PigE15B.10SOM <- FindVariableFeatures(object = PigE15B.10SOM, selection.method = "vst", nfeatures = 5000)
PigE15C.10SOM <- FindVariableFeatures(object = PigE15C.10SOM, selection.method = "vst", nfeatures = 5000)

## Scale data
PigE11.5A.EPS <- ScaleData(PigE11.5A.EPS, features = rownames(PigE11.5A.EPS))
PigE11.5B.EPS <- ScaleData(PigE11.5B.EPS, features = rownames(PigE11.5B.EPS))
PigE12A.PS <- ScaleData(PigE12A.PS, features = rownames(PigE12A.PS))
PigE12B.PS <- ScaleData(PigE12B.PS, features = rownames(PigE12B.PS))
PigE12.5A.PS2 <- ScaleData(PigE12.5A.PS2, features = rownames(PigE12.5A.PS2))
PigE12.5B.PS2 <- ScaleData(PigE12.5B.PS2, features = rownames(PigE12.5B.PS2))
PigE13A.NP <- ScaleData(PigE13A.NP, features = rownames(PigE13A.NP))
PigE13B.NP <- ScaleData(PigE13B.NP, features = rownames(PigE13B.NP))
PigE13C.NP <- ScaleData(PigE13C.NP, features = rownames(PigE13C.NP))
PigE13D.NP <- ScaleData(PigE13D.NP, features = rownames(PigE13D.NP))
PigE13E.NP <- ScaleData(PigE13E.NP, features = rownames(PigE13E.NP))
PigE13.5A.NG <- ScaleData(PigE13.5A.NG, features = rownames(PigE13.5A.NG))
PigE13.5B.NG <- ScaleData(PigE13.5B.NG, features = rownames(PigE13.5B.NG))
PigE13.5C.NG <- ScaleData(PigE13.5C.NG, features = rownames(PigE13.5C.NG))
PigE13.5D.NG <- ScaleData(PigE13.5D.NG, features = rownames(PigE13.5D.NG))
PigE14A.6SOM <- ScaleData(PigE14A.6SOM, features = rownames(PigE14A.6SOM))
PigE14B.6SOM <- ScaleData(PigE14B.6SOM, features = rownames(PigE14B.6SOM))
PigE14.5A.8SOM <- ScaleData(PigE14.5A.8SOM, features = rownames(PigE14.5A.8SOM))
PigE14.5B.8SOM <- ScaleData(PigE14.5B.8SOM, features = rownames(PigE14.5B.8SOM))
PigE14.5C.8SOM <- ScaleData(PigE14.5C.8SOM, features = rownames(PigE14.5C.8SOM))
PigE14.5D.8SOM <- ScaleData(PigE14.5D.8SOM, features = rownames(PigE14.5D.8SOM))
PigE15A.10SOM <- ScaleData(PigE15A.10SOM, features = rownames(PigE15A.10SOM))
PigE15B.10SOM <- ScaleData(PigE15B.10SOM, features = rownames(PigE15B.10SOM))
PigE15C.10SOM <- ScaleData(PigE15C.10SOM, features = rownames(PigE15C.10SOM))

#Run PCA
# We calculate 100 PCs for each sample but we will only use the ones that explain a signifficant amount of variance in the data
PigE11.5A.EPS <- RunPCA(object = PigE11.5A.EPS,  npcs = 100, verbose = FALSE)
PigE11.5B.EPS <- RunPCA(object = PigE11.5B.EPS,  npcs = 100, verbose = FALSE)
PigE12A.PS <- RunPCA(object = PigE12A.PS,  npcs = 100, verbose = FALSE)
PigE12B.PS <- RunPCA(object = PigE12B.PS,  npcs = 100, verbose = FALSE)
PigE12.5A.PS2 <- RunPCA(object = PigE12.5A.PS2,  npcs = 100, verbose = FALSE)
PigE12.5B.PS2 <- RunPCA(object = PigE12.5B.PS2,  npcs = 100, verbose = FALSE)
PigE13A.NP <- RunPCA(object = PigE13A.NP,  npcs = 100, verbose = FALSE)
PigE13B.NP <- RunPCA(object = PigE13B.NP,  npcs = 100, verbose = FALSE)
PigE13C.NP <- RunPCA(object = PigE13C.NP,  npcs = 100, verbose = FALSE)
PigE13D.NP <- RunPCA(object = PigE13D.NP,  npcs = 100, verbose = FALSE)
PigE13E.NP <- RunPCA(object = PigE13E.NP,  npcs = 100, verbose = FALSE)
PigE13.5A.NG <- RunPCA(object = PigE13.5A.NG,  npcs = 100, verbose = FALSE)
PigE13.5B.NG <- RunPCA(object = PigE13.5B.NG,  npcs = 100, verbose = FALSE)
PigE13.5C.NG <- RunPCA(object = PigE13.5C.NG,  npcs = 100, verbose = FALSE)
PigE13.5D.NG <- RunPCA(object = PigE13.5D.NG,  npcs = 100, verbose = FALSE)
PigE14A.6SOM <- RunPCA(object = PigE14A.6SOM,  npcs = 100, verbose = FALSE)
PigE14B.6SOM <- RunPCA(object = PigE14B.6SOM,  npcs = 100, verbose = FALSE)
PigE14.5A.8SOM <- RunPCA(object = PigE14.5A.8SOM,  npcs = 100, verbose = FALSE)
PigE14.5B.8SOM <- RunPCA(object = PigE14.5B.8SOM,  npcs = 100, verbose = FALSE)
PigE14.5C.8SOM <- RunPCA(object = PigE14.5C.8SOM,  npcs = 100, verbose = FALSE)
PigE14.5D.8SOM <- RunPCA(object = PigE14.5D.8SOM,  npcs = 100, verbose = FALSE)
PigE15A.10SOM <- RunPCA(object = PigE15A.10SOM,  npcs = 100, verbose = FALSE)
PigE15B.10SOM <- RunPCA(object = PigE15B.10SOM,  npcs = 100, verbose = FALSE)
PigE15C.10SOM <- RunPCA(object = PigE15C.10SOM,  npcs = 100, verbose = FALSE)


#Save all the individual objects after removal of low quality cells and doublets
SaveH5Seurat(PigE11.5A.EPS, "PigE11.5A.EPS.2", overwrite = TRUE, verbose = TRUE)
SaveH5Seurat(PigE11.5B.EPS, "PigE11.5B.EPS.2", overwrite = TRUE, verbose = TRUE)
SaveH5Seurat(PigE12A.PS, "PigE12A.PS.2", overwrite = TRUE, verbose = TRUE)
SaveH5Seurat(PigE12B.PS, "PigE12B.PS.2", overwrite = TRUE, verbose = TRUE)
SaveH5Seurat(PigE12.5A.PS2, "PigE12.5A.PS2.2", overwrite = TRUE, verbose = TRUE)
SaveH5Seurat(PigE12.5B.PS2, "PigE12.5B.PS2.2", overwrite = TRUE, verbose = TRUE)
SaveH5Seurat(PigE13A.NP, "PigE13A.NP.2", overwrite = TRUE, verbose = TRUE)
SaveH5Seurat(PigE13B.NP, "PigE13B.NP.2", overwrite = TRUE, verbose = TRUE)
SaveH5Seurat(PigE13C.NP, "PigE13C.NP.2", overwrite = TRUE, verbose = TRUE)
SaveH5Seurat(PigE13D.NP, "PigE13D.NP.2", overwrite = TRUE, verbose = TRUE)
SaveH5Seurat(PigE13E.NP, "PigE13E.NP.2", overwrite = TRUE, verbose = TRUE)
SaveH5Seurat(PigE13.5A.NG, "PigE13.5A.NG.2", overwrite = TRUE, verbose = TRUE)
SaveH5Seurat(PigE13.5B.NG, "PigE13.5B.NG.2", overwrite = TRUE, verbose = TRUE)
SaveH5Seurat(PigE13.5C.NG, "PigE13.5C.NG.2", overwrite = TRUE, verbose = TRUE)
SaveH5Seurat(PigE13.5D.NG, "PigE13.5D.NG.2", overwrite = TRUE, verbose = TRUE)
SaveH5Seurat(PigE14A.6SOM, "PigE14A.6SOM.2", overwrite = TRUE, verbose = TRUE)
SaveH5Seurat(PigE14B.6SOM, "PigE14B.6SOM.2", overwrite = TRUE, verbose = TRUE)
SaveH5Seurat(PigE14.5A.8SOM, "PigE14.5A.8SOM.2", overwrite = TRUE, verbose = TRUE)
SaveH5Seurat(PigE14.5B.8SOM, "PigE14.5B.8SOM.2", overwrite = TRUE, verbose = TRUE)
SaveH5Seurat(PigE14.5C.8SOM, "PigE14.5C.8SOM.2", overwrite = TRUE, verbose = TRUE)
SaveH5Seurat(PigE14.5D.8SOM, "PigE14.5D.8SOM.2", overwrite = TRUE, verbose = TRUE)
SaveH5Seurat(PigE15A.10SOM, "PigE15A.10SOM.2", overwrite = TRUE, verbose = TRUE)
SaveH5Seurat(PigE15B.10SOM, "PigE15B.10SOM.2", overwrite = TRUE, verbose = TRUE)
SaveH5Seurat(PigE15C.10SOM, "PigE15C.10SOM.2", overwrite = TRUE, verbose = TRUE)

# Add load option as needed
PigE11.5A.EPS <- LoadH5Seurat("/Users/lukesimpson/Desktop/Integrated analysis/PigE11.5A.EPS.h5seurat", reductions = TRUE,meta.data = TRUE)
PigE11.5B.EPS <- LoadH5Seurat("/Users/lukesimpson/Desktop/Integrated analysis/PigE11.5B.EPS.h5seurat", reductions = TRUE,meta.data = TRUE)
PigE12A.PS <- LoadH5Seurat("/Users/lukesimpson/Desktop/Integrated analysis/PigE12A.PS.h5seurat", reductions = TRUE,meta.data = TRUE)
PigE12B.PS <- LoadH5Seurat("/Users/lukesimpson/Desktop/Integrated analysis/PigE12B.PS.h5seurat", reductions = TRUE,meta.data = TRUE)
PigE12.5A.PS2 <- LoadH5Seurat("/Users/lukesimpson/Desktop/Integrated analysis/PigE12.5A.PS2.h5seurat", reductions = TRUE,meta.data = TRUE)
PigE12.5B.PS2 <- LoadH5Seurat("/Users/lukesimpson/Desktop/Integrated analysis/PigE12.5B.PS2.h5seurat", reductions = TRUE,meta.data = TRUE)
PigE13A.NP <- LoadH5Seurat("/Users/lukesimpson/Desktop/Integrated analysis/PigE13A.NP.h5seurat", reductions = TRUE,meta.data = TRUE)
PigE13B.NP <- LoadH5Seurat("/Users/lukesimpson/Desktop/Integrated analysis/PigE13B.NP.h5seurat", reductions = TRUE,meta.data = TRUE)
PigE13C.NP <- LoadH5Seurat("/Users/lukesimpson/Desktop/Integrated analysis/PigE13C.NP.h5seurat", reductions = TRUE,meta.data = TRUE)
PigE13D.NP <- LoadH5Seurat("/Users/lukesimpson/Desktop/Integrated analysis/PigE13D.NP.h5seurat", reductions = TRUE,meta.data = TRUE)
PigE13E.NP <- LoadH5Seurat("/Users/lukesimpson/Desktop/Integrated analysis/PigE13E.NP.h5seurat", reductions = TRUE,meta.data = TRUE)
PigE13.5A.NG <- LoadH5Seurat("/Users/lukesimpson/Desktop/Integrated analysis/PigE13.5A.NG.h5seurat", reductions = TRUE,meta.data = TRUE)
PigE13.5B.NG <- LoadH5Seurat("/Users/lukesimpson/Desktop/Integrated analysis/PigE13.5B.NG.h5seurat", reductions = TRUE,meta.data = TRUE)
PigE13.5C.NG <- LoadH5Seurat("/Users/lukesimpson/Desktop/Integrated analysis/PigE13.5C.NG.h5seurat", reductions = TRUE,meta.data = TRUE)
PigE13.5D.NG <- LoadH5Seurat("/Users/lukesimpson/Desktop/Integrated analysis/PigE13.5D.NG.h5seurat", reductions = TRUE,meta.data = TRUE)
PigE14A.6SOM <- LoadH5Seurat("/Users/lukesimpson/Desktop/Integrated analysis/PigE14A.6SOM.h5seurat", reductions = TRUE,meta.data = TRUE)
PigE14B.6SOM <- LoadH5Seurat("/Users/lukesimpson/Desktop/Integrated analysis/PigE14B.6SOM.h5seurat", reductions = TRUE,meta.data = TRUE)
PigE14.5A.8SOM <- LoadH5Seurat("/Users/lukesimpson/Desktop/Integrated analysis/PigE14.5A.8SOM.h5seurat", reductions = TRUE,meta.data = TRUE)
PigE14.5B.8SOM <- LoadH5Seurat("/Users/lukesimpson/Desktop/Integrated analysis/PigE14.5B.8SOM.h5seurat", reductions = TRUE,meta.data = TRUE)
PigE14.5C.8SOM <- LoadH5Seurat("/Users/lukesimpson/Desktop/Integrated analysis/PigE14.5C.8SOM.h5seurat", reductions = TRUE,meta.data = TRUE)
PigE14.5D.8SOM <- LoadH5Seurat("/Users/lukesimpson/Desktop/Integrated analysis/PigE14.5D.8SOM.h5seurat", reductions = TRUE,meta.data = TRUE)
PigE15A.10SOM <- LoadH5Seurat("/Users/lukesimpson/Desktop/Integrated analysis/PigE15A.10SOM.h5seurat", reductions = TRUE,meta.data = TRUE)
PigE15B.10SOM <- LoadH5Seurat("/Users/lukesimpson/Desktop/Integrated analysis/PigE15B.10SOM.h5seurat", reductions = TRUE,meta.data = TRUE)
PigE15C.10SOM <- LoadH5Seurat("/Users/lukesimpson/Desktop/Integrated analysis/PigE15C.10SOM.h5seurat", reductions = TRUE,meta.data = TRUE)

# select features that are repeatedly variable across datasets for integration
featuresAll <- SelectIntegrationFeatures(object.list = c(PigE11.5A.EPS,	PigE11.5B.EPS,	PigE12A.PS,	PigE12B.PS,	PigE12.5A.PS2,	PigE12.5B.PS2,	PigE13A.NP,	PigE13B.NP,	PigE13C.NP,	PigE13D.NP,	PigE13E.NP,	PigE13.5A.NG,	PigE13.5B.NG,	PigE13.5C.NG,	PigE13.5D.NG,	PigE14A.6SOM,	PigE14B.6SOM,	PigE14.5B.8SOM,	PigE14.5C.8SOM,	PigE14.5D.8SOM,	PigE15A.10SOM,	PigE15B.10SOM,	PigE15C.10SOM), nfeatures = 4000) #was 2000
# select features for downstream integration, and run PCA on each object in the list, which is required for running the alternative reciprocal PCA workflow
Pig.anchorsAll <- FindIntegrationAnchors(object.list = c(PigE11.5A.EPS,	PigE11.5B.EPS,	PigE12A.PS,	PigE12B.PS,	PigE12.5A.PS2,	PigE12.5B.PS2,	PigE13A.NP,	PigE13B.NP,	PigE13C.NP,	PigE13D.NP,	PigE13E.NP,	PigE13.5A.NG,	PigE13.5B.NG,	PigE13.5C.NG,	PigE13.5D.NG,	PigE14A.6SOM,	PigE14B.6SOM,	PigE14.5B.8SOM,	PigE14.5C.8SOM,	PigE14.5D.8SOM,	PigE15A.10SOM,	PigE15B.10SOM,	PigE15C.10SOM), reference = c(1,4,5,10,12,16,18,21), anchor.features = featuresAll, reduction = "rpca", dims = 1:50) #
# Perform weighted integration
Pig.all.combined <- IntegrateData(anchorset = Pig.anchorsAll, new.assay.name = "integrated")
#cell cycle scoring 
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
Pig.all.combined <- CellCycleScoring(Pig.all.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# Rescale the integrated data
Pig.all.combined <- ScaleData(Pig.all.combined, features = rownames(Pig.all.combined), verbose = TRUE, )
Pig.all.combined <- RunPCA(Pig.all.combined, npcs = 100, verbose = FALSE)
Pig.all.combined <- FindNeighbors(Pig.all.combined, reduction = "pca", dims = 1:50)
# Cluster resolution was chosen arbitrarily based on the resolution that gives a number closest to the number of cell types we are expecting
Pig.all.combined <- FindClusters(Pig.all.combined, resolution = 2, group.singletons = TRUE)
Pig.all.combined <- RunUMAP(Pig.all.combined, reduction = "pca", dims = 1:50, n.neighbors = 50L, seed.use = 1, min.dist = 0.4, return.model = TRUE) #spread = 1, n.components = 3L,
# Check for high MT clusters that span multiple cell types or low count clusters
DimPlot(Pig.all.combined, reduction ="umap", group.by= "seurat_clusters",label = TRUE) + NoLegend()
FeaturePlot(Pig.all.combined, features = "nFeature_RNA", reduction = "umap", order = T)
# Remove poor quality clusters (Clusters which span across multiple lineages and have high mito expression or clusters with low counts)
Pig.all.combined <- subset(Pig.all.combined, idents = c("31","33"), invert = TRUE)
Pig.all.combined <- subset(Pig.all.combined, idents = c("5"), invert = TRUE)

#Re-run PCA, umap and NN finding and clustering
Pig.all.combined <- ScaleData(Pig.all.combined, verbose = TRUE)
Pig.all.combined <- RunPCA(object = Pig.all.combined,  npcs = 30, verbose = FALSE, features = NULL) # features = NULL and rownames(Pig.int) are similar and decent
Pig.all.combined <- FindNeighbors(Pig.all.combined, reduction = "pca", dims = 1:25)
Pig.all.combined <- FindClusters(Pig.all.combined, resolution = 1.6, group.singletons = TRUE) #1.6
Pig.all.combined <- RunUMAP(Pig.all.combined, reduction = "pca", dims = 1:25, n.neighbors = 30L, min.dist = 0.5, return.model = TRUE,  seed.use = 2, umap.method = "umap-learn", metric = "correlation") # 1:20/21 is best, 22 is ok, 23 is same, 1:19 great,  n.epochs = 200 is good at 0.5 dist umap.method = "umap-learn", metric = "correlation" , dims 1:22 seed 42 best, min dist 0.6 good 1 PGCs are in good place, umap.method = "umap-learn", metric = "correlation" #spread = 1, n.components = 3L,

# Save integrated seurat objects
SaveH5Seurat(Pig.all.combined, "Pig.all.combined", overwrite = TRUE, verbose = TRUE)
#Load 
Pig.all.combined <- LoadH5Seurat("/Users/lukesimpson/Pig.all.combined.h5seurat")

