## Annotation of cell types for single-cell atlas of pig gastrulation ##
# Load libraries
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
  library(viridis)  # For colour-blind friendly plots
}

Pig.all.combined <- LoadH5Seurat("/Users/lukesimpson/Pig.all.combined.h5seurat")

DefaultAssay(Pig.all.combined) <- "integrated"
# Does this replicate the clustering last thing tested
Pig.all.combined <- ScaleData(Pig.all.combined, verbose = TRUE)
Pig.all.combined <- RunPCA(object = Pig.all.combined,  npcs = 30, verbose = FALSE, features = NULL) 
Pig.all.combined <- FindNeighbors(Pig.all.combined, reduction = "pca", dims = 1:25)
Pig.all.combined <- FindClusters(Pig.all.combined, resolution = 1.6, group.singletons = TRUE) #1.6 gives 35 clusters
Pig.all.combined <- RunUMAP(Pig.all.combined, reduction = "pca", dims = 1:25, n.neighbors = 30L, min.dist = 0.5, return.model = TRUE,  seed.use = 2, umap.method = "umap-learn", metric = "correlation") 

#Query marker genes and decide on cluster ID 
#The new seurat object is too big for query mapping or integration on a personal computer
Markers <- FindAllMarkers(Pig.all.combined)
write.csv(Markers, "Markers.csv")

# Plot a heat map of the top 10 markers for each cluster NOTE: This didn't appear in the final paper
Top_10_markers <- Markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 10) %>% ungroup() 
# Plot a heat map
DoHeatmap(Pig.all.combined, features = Top_10_markers$gene) + NoLegend() 

# Lets look at the expression of the marker genes used by Pijuan-Sala et al., 2019 NOTE:This didn't appear in the final paper
Mouse_Diagnostic_genes <- c("ENSSSCG00000001393",	"DNMT3B",	"EPCAM",	"NANOG",	"UTF1",	"EOMES",	"NKX1-2",	"IFITM3",	"DND1",	"DPPAЗ",	"HHEX",	"OTX2",	"MIXL1",	"GSC",	"LHX1",	"FGF8",	"PAX7",	"FOXA2",	"CHRD",	"NOTO",	"CER1",	"APELA",	"MESP1",	"LEFTY2",	"MESP2",	"OSR1",	"CDX2",	"HES7",	"HOXB1",	"CDX1",	"GBX2",	"TBX1",	"MEOX1",	"TCF15",	"ALDH1A2",	"DLL1",	"TBX6",	"ISL1",	"TCF21",	"SMARCD3",	"ACTA2",	"TAGLN",	"MYL7",	"NKX2-5",	"MYL4",	"TNNT2",	"HOXA10",	"HOXA11",	"TBX4",	"BMP4",	"AHNAK",	"PMP22",	"KRT18",	"KRT8",	"IGF2",	"ETV2",	"KDR",	"ANXA5",	"PECAM1",	"CDH5",	"LMO2",	"RUNX1",	"CITED4",	"GATA1",	"HBBG1",	"GYPA",	"HBA1",	"HBA2",	"HES3",	"EPHA5",	"CDX4",	"HOXB9",	"SOX2",	"IRX3",	"HESX1",	"SIX3",	"SOX9",	"PAX3",	"TFAP2A",	"FOXD3",	"DLX2",	"SOX10",	"EN1",	"PAX2",	"PAX6",	"FOXG1",	"TRP63",	"GRHL2",	"GRHL3",	"DKK1",	"KRT19",	"AMOT",	"SPINK1",	"EMB",	"CYSTM1",	"APOE",	"АРОА2",	"TTR",	"TFAP2C",	"ASCL2",	"ELF5",	"SPARC",	"PLAT",	"LAMB1")

Smallavgexp <- AverageExpression(Pig.all.combined, return.seurat = T)
DoHeatmap(object = Smallavgexp, assay = "RNA", draw.lines = F, size = 2, disp.min = -0.5, group.by = "ident", features = as.character(Mouse_Diagnostic_genes)) + scale_fill_gradientn(colors = viridis(50))+ theme(text = element_text(size = 10)) + theme(legend.position = "none")

# Now lets plot a heat map showing the expression of known cell type markers in each cluster 
Celltype.Diagnostic.markers <- c("NANOG", "PRDM1", "DND1", "ENSSSCG00000001393", "DNMT3B", "EPCAM", "FST","EOMES","CDX2", "TBXT", "LHX1", "CER1", "GSC", "CHRD", "DKK1", "FOXA2", "SOX17", "LEFTY2", "MESP1", "HOXB1", "TBX6", "DLL1", "MEOX1", "TCF15", "FOXC1", "FOXC2", "TBX1", "MYF5", "OSR1", "HAND1", "TCF21", "ISL1", "SFRP5", "ACTA2", "TAGLN", "MYL4", "HOXA11", "HOXA10", "BMP4", "POSTN", "TMEM88", "COL3A1", "GATA6", "ETV2", "KDR", "CDH5", "LMO2", "RUNX1", "PECAM1", "GATA1", "PAX6", "PAX7", "TFAP2C", "GRHL3", "TFAP2A","HOXB9", "GATA3","GABRP","HEY1","DAB2", "AMOT", "GATA4", "APOE", "VIL1", "TTR", "HNF1B")
DoHeatmap(object = Smallavgexp, assay = "RNA", draw.lines = F, size = 2, disp.min = -0.5, group.by = "ident", features = as.character(Celltype.Diagnostic.markers)) + scale_fill_gradientn(colors = viridis(50))+ theme(text = element_text(size = 10)) + theme(legend.position = "none")

# We can also use HOX genes as a sanity check (later cell types should have a HOX profile consistent with their identity) NOTE: This didn't appear in the final paper
Anterior_posterior_genes <- c("HOXA1","HOXB1", "HOXD1","HOXA2","HOXB2","HOXA3","HOXB3","HOXD3","HOXB4","HOXC4","HOXA5","HOXB5","HOXC5","HOXB6","HOXC6","HOXA7","HOXB7","HOXB8","HOXC8","HOXD8","HOXB9","HOXC9","HOXD9","HOXA10","HOXD10","HOXA11","HOXD11","HOXC12","HOXD12","HOXA13","HOXB13","HOXD13")
DoHeatmap(object = Smallavgexp, assay = "RNA", draw.lines = F, size = 2, disp.min = -0.5, group.by = "ident", features = as.character(Anterior_posterior_genes)) + scale_fill_gradientn(colors = viridis(50))+ theme(text = element_text(size = 10)) + theme(legend.position = "none")

# We can now name our clusters
Pig.all.combined <- SetIdent(Pig.all.combined, value = Pig.all.combined@meta.data[["seurat_clusters"]])
Pig.all.combined <- RenameIdents(Pig.all.combined, 
                    "0"=                     "Nascent Mesoderm 1",
                    "1"=                     "Mesenchyme 3",
                    "2"=                     "Somitic Mesoderm",
                    "3"=                     "Posterior Mixed mesoderm",
                    "4"=                    "Epiblast 1",
                    "5"=                     "Intermediate Mesoderm 2",
                    "6"=                     "ExE Mesoderm 1",
                    "7"=                     "Primitive streak 1",
                    "8"=                     "ExE Endoderm",
                    "9"=                     "ExE Mesoderm 2",
                    "10"=                     "Pharyngeal Mesoderm 1",
                    "11"=                     "Nascent Mesoderm 2",
                    "12"=                     "Amnion",
                    "13"=                     "Primitive Streak 2",
                    "14"=                     "Intermediate Mesoderm 1",
                    "15"=                     "Mixed Mesoderm",
                    "16"=                     "Anterior Surface Ectoderm",
                    "17"=                     "Definitive Endoderm",
                    "18"=                     "Epiblast 4",
                    "19"=                     "Mesenchyme 2",
                    "20"=                     "Mesenchyme 1",
                    "21"=                     "Posterior Surface Ectoderm",
                    "22"=                     "Gut/Hypoblast",
                    "23"=                     "Epiblast 2",
                    "24"=                     "Primitive Streak 3",
                    "25"=                     "Epiblast 3",
                    "26"=                     "Brain",
                    "27"=                     "Anterior Primitive Streak/Node",
                    "28"=                     "Haematoendothelial Progenitors",
                    "29"=                     "Trophoblast",
                    "30"=                     "Allantois", 
                    "31"=                     "Spinal Cord",
                    "32"=                     "Presomitic Mesoderm",
                    "33"=                     "Blood Progenitors",
                    "34"=                     "Cardiac Mesoderm",
                    "35"=                     "PGC")

Pig.all.combined$CellType <- Idents(Pig.all.combined)

# Now we replot the heatmap with our new cell identities and custom colours
celltype_colors <- c("Trophoblast"=	"#3565A7","Lateral Plate Mesoderm"= "#E0057A", "Intermediate mesoderm 2"= "#f4192b", "Amnion"=	"#5891BF","Anterior Surface Ectoderm"=	"#D3F1FD", "Posterior Surface Ectoderm"="#C1DDBB","Caudal Epiblast"=	"#92CBC1","Spinal Cord" = "#2B645E","Epiblast 4"=	"#CAD2CE",  "Epiblast 1"=	"#D5BF9E","Epiblast 2"=	"#E8D6A8","Epiblast 3"=	"#F8EAB0","Anterior Primitive Streak/Node"=	"#FFF478","Definitive Endoderm"=	"#E6FFA6", "Gut/Hypoblast"=	"#AFE76A", "ExE Endoderm"=	"#327D21","Primitive streak 1"=	"#FFE173","Primitive Streak 3"=	"#FFCD8B", "Primitive Streak 2"=	"#FBAD54", "Nascent Mesoderm 2"=	"#F99878","Nascent Mesoderm 1"=	"#EE7C54", "Intermediate Mesoderm 1"= "#ff624a", "Posterior Mixed mesoderm"=	"#B81828","ExE Mesoderm 2"=	"#A2545E", "ExE Mesoderm 1"=	"#821825", "Mesenchyme 3"=	"#724545","Mesenchyme 2"=	"#430707","Cardiac Mesoderm"=	"#000000", "Allantois" = "#AE1667","Pharyngeal Mesoderm"=	"#ee71b5","Presomitic Mesoderm"=	"#E076C5", "Somitic Mesoderm"=	"#F9D1EE","Haematoendothelial Progenitors"=	"#8F4986","Blood Progenitors"=	"#590D57","Mesenchyme 1"=	"#1b0840", "PGC"=	"#8249B4")
Pig.all.combined <- SetIdent(Pig.all.combined, value = Pig.all.combined@meta.data[["CellType"]])
Smallavgexp <- AverageExpression(Pig.all.combined, return.seurat = T)
DoHeatmap(object = Smallavgexp, assay = "RNA", draw.lines = F, size = 2, disp.min = -0.5, group.by = "ident", group.colors =celltype_colors, features = as.character(Celltype.Diagnostic.markers)) + scale_fill_gradientn(colors = viridis(50))+ theme(text = element_text(size = 10)) + theme(legend.position = "none")

# Extract the data for the source data file
average_expression_matrix <- as.data.frame(Smallavgexp[["RNA"]]@scale.data)
average_expression_matrix <- average_expression_matrix[Celltype.Diagnostic.markers,]
average_expression_matrix$Gene <- rownames(average_expression_matrix)
# Write the data to a CSV file
write.csv(average_expression_matrix, "C:/Users/sbzlas1/Documents/Luke Simpson Remote computer  temporary files/SupplementaryFigure1e.csv", row.names = FALSE)

# Save annotated pig atlas
SaveH5Seurat(Pig.all.combined, "Pig.all.combined", overwrite = TRUE, verbose = TRUE)
#Load 
Pig.all.combined <- LoadH5Seurat("/Users/lukesimpson/Pig.all.combined.h5seurat")
