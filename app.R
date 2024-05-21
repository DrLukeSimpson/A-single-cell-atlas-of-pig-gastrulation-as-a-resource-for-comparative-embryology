### Alberio lab shiny app ###
#install.packages("shiny")
#install.packages("recolorize")
#install.packages('rsconnect')

#rsconnect::setAccountInfo(name='luke-simpson-university-of-nottingham',token='2054787240274AD2A047053BFA961A7E', secret='zifUP1MNMx42mzDV4MfqYmuwwso0dUthY0NrGMXF')

#options(rsconnect.max.bundle.size=3772417930)

library(shiny)
library(Seurat)
library(dplyr)
library(ggplot2)
library(viridis)
library(tidyr)
library(magick)
library(rsconnect)

#rsconnect::deployApp('/Users/lukesimpson/Desktop/Pig developmental resource shiny app/App folder/')

#Set wd
#setwd("/Users/lukesimpson/Desktop/Pig developmental resource shiny app/App folder/")
# Load objects
seurat_object <- readRDS(file = "Piglargeseuratobjectforshinyappnointegrated09042024.RDS")
seurat_object2 <- readRDS(file = "E5-E10Seurat.RDS")

# Ensure RNA is default assay for plotting
DefaultAssay(seurat_object) <- "RNA"
DefaultAssay(seurat_object2) <- "RNA"
celltype_colors <- c("Trophoblast"=	"#3565A7","Lateral Plate Mesoderm"= "#E0057A", "Intermediate mesoderm 2"= "#f4192b", "Amnion"=	"#5891BF","Anterior Surface Ectoderm"=	"#D3F1FD", "Posterior Surface Ectoderm"="#C1DDBB","Caudal Epiblast"=	"#92CBC1","Spinal Cord" = "#2B645E","Epiblast 4"=	"#CAD2CE",  "Epiblast 1"=	"#D5BF9E","Epiblast 2"=	"#E8D6A8","Epiblast 3"=	"#F8EAB0","Anterior Primitive Streak/Node"=	"#FFF478","Definitive Endoderm"=	"#E6FFA6", "Gut/Hypoblast"=	"#AFE76A", "ExE Endoderm"=	"#327D21","Primitive streak 1"=	"#FFE173","Primitive Streak 3"=	"#FFCD8B", "Primitive Streak 2"=	"#FBAD54", "Nascent Mesoderm 2"=	"#F99878","Nascent Mesoderm 1"=	"#EE7C54", "Intermediate Mesoderm 1"= "#ff624a", "Posterior Mixed mesoderm"=	"#B81828","ExE Mesoderm 2"=	"#A2545E", "ExE Mesoderm 1"=	"#821825", "Mesenchyme 3"=	"#724545","Mesenchyme 2"=	"#430707","Cardiac Mesoderm"=	"#000000", "Allantois" = "#AE1667","Pharyngeal Mesoderm"=	"#ee71b5","Presomitic Mesoderm"=	"#E076C5", "Somitic Mesoderm"=	"#F9D1EE","Haematoendothelial Progenitors"=	"#8F4986","Blood Progenitors"=	"#590D57","Mesenchyme 1"=	"#1b0840", "PGC"=	"#8249B4")

# Convert the colors to HSV and extract the hue
hues <- sapply(celltype_colors, function(hex) {
  rgb <- col2rgb(hex) / 255
  hsv <- rgb2hsv(rgb[1,], rgb[2,], rgb[3,])
  return(hsv[1,]) # Extract the hue component
})

# Order the celltype_colors by hue
# Filter out NAs just in case they occur
valid_indices <- !is.na(hues)
ordered_colors <- celltype_colors[valid_indices][order(hues[valid_indices])]

cellsubtype_colors <- c("ExE Endoderm"="grey99","Definitive Endoderm"="white","Epiblast 4"="#94c9bd","Epiblast 5"= "#6e6352","Early Caudal Epiblast 1"=	"#bd9c39","Early Caudal Epiblast 2" ="#78621f", "Anterior Primitive Streak"="#FFF478","Epiblast 1"="#D5BF9E", "Nascent Mesoderm 1"="#EE7C54", "Nascent Mesoderm 2"="#F99878", "Nascent Mesoderm 3"="#F4192B", "Nascent Mesoderm 4"="#f75f2d", "Node"="#d6d672", "Primitive Streak 1"="#FFE173", "Primitive Streak 2"="#FBAD54", "Primitive Streak 3"="#FFCD8B","Allantois"="#AE1667", "Amnion"="#5891BF", "Anterior Hypoblast/AVE"="#45b069", "Anterior Somitic Mesoderm 1"="#f9d1ee", "Anterior Somitic Mesoderm 2"="#e076c5", "Anterior Somitic Mesoderm 3"="#c892be", "Anterior Surface Ectoderm"="#D3F1FD", "Blood Progenitors"="#590D57", "Cardiac Mesoderm"="#000000", "Caudal Epiblast"="lightgrey", "Cranial Mesoderm"="#1b0845", "Definitive Endoderm/Foregut"="#E6FFA6", "Definitive Endoderm/Hindgut"="#bdff59", "Dermomyotome/Sclerotome"="#5D1A60", "Epiblast 2"="#E8D6A8", "Epiblast 3"="#F8EAB0", "ExE APOP"="#55695b", "ExE Endo 1"="#327D35", "ExE Endo 2"="#327D21", "ExE Endo 3"="#122e0b", "ExE Endo 4"="#9bc7a8", "ExE Endo 5"="#4bd444", "ExE Endo 6"="#82ffdc", "ExE Mesoderm 1"="#821825", "ExE Mesoderm 2"="#A2545E", "Haematoendothelial Progenitors"="#8F4970", "Hypoblast 1"="#69cf77", "Hypoblast 2"="#438f4d", "InterHypo"="#6dd955", "Intermediate Mesoderm 1"="#ff624a", "Intermediate mesoderm 2"="#EE71B6", "Lateral Plate Mesoderm"="#cc3300", "Mesenchyme 1"="#1b0840", "Mesenchyme 2"="#430707", "Mesenchyme 3"="#724545", "Midgut 1"="#8ae34b", "Midgut 2"="#8aa83d","PGC"="#8249B4", "Pharyngeal Mesoderm"="#E0057A", "Posterior Mixed mesoderm"="#B81828", "Posterior Somitic Mesoderm"="#e0057a", "Posterior Surface Ectoderm"="#C1DDBB", "Presomitic Mesoderm"="#8F4986", "Spinal Cord"="#2B645E", "Trophoblast"="#3565A7", "YS Endo"="#55d98e")

hues2 <- sapply(cellsubtype_colors, function(hex) {
  rgb <- col2rgb(hex) / 255
  hsv <- rgb2hsv(rgb[1,], rgb[2,], rgb[3,])
  return(hsv[1,]) # Extract the hue component
})

valid_indices2 <- !is.na(hues2)
ordered_colors2 <- cellsubtype_colors[valid_indices2][order(hues2[valid_indices2])]
ordered_colors3 <- c("Hypoblast_e11"="#f4d839","Morula_e5"="#dc0100", "Epiblast_e11"="#e87b54","ICM_e6"="#5975fa","Trophectoderm_e6"="#a0cde9","Hypoblast_e8"="#a9cc45","Epiblast_e8"="#b1ee95", "Morula_e8"="#dc0100", "Trophectoderm_e8"="#a0cde9")

# Assume 'seurat_object' is your loaded Seurat object
ui <- navbarPage(
  tags$head(
    tags$style(HTML("
      @font-face {
        font-family: 'Acumin Pro';
        src: url('www/acumin-pro/Acumin-RPro.otf') format('opentype');
        font-weight: normal;
        font-style: normal;
      }
      @font-face {
        font-family: 'Acumin Pro';
        src: url('www/acumin-pro/Acumin-BdPro.otf') format('opentype');
        font-weight: bold;
        font-style: normal;
      }
      .navbar-default .navbar-brand {
        font-size: 30px;
        color: black !important;
        font-family: 'Acumin Pro', sans-serif;
        font-weight: bold;
      }
      .footer-text {
        position: fixed;
        width: 100%;
        text-align: center;
        bottom: 0;
        color: grey65;
        background-color: white; /* Ensures visibility over any page content */
        padding: 10px 0;
      }
    "))
  ),
  title = "Pig Developmental Resource",
  
  tabPanel("Landing Page",
           fluidPage(
             fluidRow(
               column(6,
                      tags$h2("Welcome"),
                      tags$div(style="font-size: 20px; font-weight: bold;", 
                               "to the Pig Developmental Resource, here you will  find resources relating to the study of pig embryology including the pig gastrulation atlas."),
                      tags$br(),  # line break
                      tags$p(tags$b("Gastulation dataset overview:"), "UMAP plots of the data can be viewed highlighting cell types or cell subtypes as they appear in Simpson et al., 2023. The whole dataset covering E11.5 to E15 or individual stages can be viewed."),
                      tags$p(tags$b("Gene Expression:"), "Normalised counts of gene expression across E11.5 to E15 pig embryos from Simpson et al., 2023."),
                      tags$p(tags$b("Cell-type Markers:"), "Genes that demarcate a particular cell type are shown here."),
                      tags$p(tags$b("Cell Subtype Markers:"), "Genes that demarcate a particular cell type are shown here."),
                      tags$p(tags$b("Cell-type Discrimination:"), "Genes that demarcate a particular cell type are shown here."),
                      tags$p(tags$b("Pre-gastrula Gene Expression:"), "Normalised counts of gene expression across E5 to E10 pig embryos from Ramos-Ibeas et al., 2019."),
                      tags$p(tags$b("Pig Resources:"), "Highlighted papers from the pig embryology community, a great starting place for anyone looking to become acquainted with this area of research.")
               ),
               column(6, img(src = "PIGSCART2.png", height = "100%", width = "100%"))
             ),
             tags$footer(tags$div(class = "footer-text", "To report any issues please contact Luke Simpson at Luke.Simpson2@nottingham.ac.uk"))
           )
  ),
  tabPanel("Dataset Overview",
           fluidPage(
             selectInput("projectionType", "Projection type", choices = c("UMAP", "tSNE")),
             selectInput("clustering", "Clustering", choices = c("Cell types", "Cell subtypes")),
             selectInput("cellSubset", "Cell subset", choices = unique(seurat_object$Stage)),
             plotOutput("dimPlot", width = "1200px", height = "450px"), 
             htmlOutput("pdfDisplay"),
             tags$footer(
               tags$div(class = "footer-text", "Please cite Simpson et al., 2023")
           ))),
  tabPanel("Gene Expression (E11.5-15)",
           fluidPage(
             selectInput("geneSelection", "Gene", choices = rownames(seurat_object@assays$RNA@data)),
             selectInput("plotType", "Plot Type", choices = c("Feature Plot", "Violin Plot", "Box and whisker plot")),
             fluidRow(
               column(6, plotOutput("genePlot", width = "490px",height = "450px")),
               column(6, plotOutput("cellTypeDimPlot", width = "1095px", height = "450px"), 
                      # Add a download button for saving the plot
                      downloadButton("downloadPlotBtn", "Download Plot"),
                      # Add a title to the sidebar
                      tags$h3("Download options"),
                      # Add a radio button for selecting the file format
                      radioButtons("fileFormat", "Select File Format:",
                                   choices = c("PDF", "PNG"),
                                   selected = "PNG"),
               tags$footer(
                 tags$div(class = "footer-text", "Please cite Simpson et al., 2023")
             )
           )))),
  tabPanel("Cell-type Markers",
           fluidPage(
             selectInput("celltype", "Cell type", choices = unique(seurat_object$Celltype)),
             plotOutput("Findmarkers1", width = "850px", height = "400px"), 
             tags$footer(
               tags$div(class = "footer-text", "Please cite Simpson et al., 2023")
           ))),
  tabPanel("Cell Subtype Markers",
           fluidPage(
             selectInput("cellsubtype", "Cell Subtype", choices = unique(seurat_object$Cellsubclusters)),
             plotOutput("Findmarkers2", width = "850px", height = "400px"), 
             tags$footer(
               tags$div(class = "footer-text", "Please cite Simpson et al., 2023")
           ))),
  tabPanel("Cell-type Discrimination",
           fluidPage(
             selectInput("celltype1", "Cell type 1", choices = unique(seurat_object$Celltype)),
             selectInput("celltype2", "Cell type 2", choices = unique(seurat_object$Celltype)),
             plotOutput("Findmarkers3", width = "850px", height = "400px"), 
             tags$footer(
               tags$div(class = "footer-text", "Please cite Simpson et al., 2023")
           ))),
  tabPanel("Gene Expression (E4-10)",
           fluidPage(
             selectInput("geneSelection2", "Gene", choices = rownames(seurat_object2@assays[["RNA"]])),
             selectInput("plotType2", "Plot Type", choices = c("Feature Plot", "Violin Plot", "Box and whisker plot")),
             fluidRow(
               column(6, plotOutput("genePlot2", width = "490px",height = "450px")),
               column(6, plotOutput("cellTypeDimPlot2", width = "835px", height = "450px")),
               tags$footer(
                 tags$div(class = "footer-text", "Please cite Ramos-Ibeas et al., 2019"))
           ))),
           
  tabPanel("Pig Resources",
           fluidPage(
             fluidRow(
               column(6,
                      tags$p(tags$b("Pluripotency and X chromosome dynamics revealed in pig pre-gastrulating embryos by single cell analysis"), "https://www.nature.com/articles/s41467-019-08387-8"),
                      tags$p(tags$b("Principles of early human development and germ cell program from conserved model systems"), "https://www.nature.com/articles/nature22812"),
                      tags$p(tags$b("Generation and characterization of stable pig pregastrulation epiblast stem cell lines"), "https://www.nature.com/articles/s41422-021-00592-9"),
                      tags$p(tags$b("A single-cell atlas of pig gastrulation as a resource for comparative embryology"), "https://www.biorxiv.org/content/10.1101/2023.08.31.555712v1"),
                      tags$p(tags$b("Cell specification and functional interactions in the pig blastocyst inferred from single-cell transcriptomics and uterine fluids proteomics"), "https://www.sciencedirect.com/science/article/pii/S0888754323002240?via%3Dihub"),
                      tags$p(tags$b("Specification and epigenomic resetting of the pig germline exhibit conservation with the human lineage"), "https://www.cell.com/cell-reports/fulltext/S2211-1247(21)00048-6?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2211124721000486%3Fshowall%3Dtrue"),
                      tags$p(tags$b("Conserved and divergent expression patterns of markers of axial development in eutherian mammals"), "https://anatomypubs.onlinelibrary.wiley.com/doi/10.1002/dvdy.24352"),
                      tags$p(tags$b("Paracrine effects of embryo-derived FGF4 and BMP4 during pig trophoblast elongation"), "https://www.sciencedirect.com/science/article/pii/S0012160614000232?via%3Dihub"),
                      tags$p(tags$b("Three-Dimensional ImmunohistochemicalCharacterization of Lineage Commitment byLocalization of T and FOXA2 in PorcinePeri-implantation Embryos"), "https://anatomypubs.onlinelibrary.wiley.com/doi/10.1002/dvdy.22602"))
                      )
             ))
 
)


server <- function(input, output) {
  
  output$dimPlot <- renderPlot({
    # Initialize variables for legend customization
    key_size <- unit(0.75, "cm")  # Key size
    key_width <- unit(0.5, "cm")  # Key width
    legend_background <- "white"  # Legend background color
    legend_key_fill <- "white"  # Legend key fill
    legend_spacing_y <- unit(0.01, "cm")  # Adjust the vertical space between legend keys
    legend_margin <- margin(t = 10, r = 295, b = 10, l = 10, unit = "pt")
    text_size <- 14  # Default text size
    
    # Determine which grouping and colors to use based on user input
    if (input$clustering == "Cell types") {
      grouping_variable <- "Celltype"
      colors_to_use <- ordered_colors
    } else if (input$clustering == "Cell subtypes") {
      grouping_variable <- "Cellsubclusters"
      colors_to_use <- ordered_colors2
      key_size <- unit(0.75, "cm")  # Adjust key size for more groups
      key_width <- unit(0.5, "cm")  # Maintain wider keys
      legend_margin <- margin(t = 10, r = 0, b = 10, l = 0, unit = "pt")
      text_size <- 12
    }
    
    # Ensure the grouping variable is present in the Seurat object
    if (!is.null(grouping_variable) && !is.null(colors_to_use)) {
      # Plot with the specified settings
      p <- DimPlot(
        object = seurat_object,
        reduction = "umap",
        group.by = grouping_variable,
        cols = colors_to_use,  # Apply the selected color palette
        pt.size = 0.75
      ) + theme(
        legend.position = "right",
        legend.text = element_text(size = text_size),
        legend.background = element_rect(fill = legend_background),  # Custom legend background
        legend.key = element_rect(fill = legend_key_fill, color = NA),  # Custom legend key appearance
        legend.key.size = key_size,  # Custom legend key size
        legend.key.width = key_width,  # Custom legend key width
        legend.spacing.y = legend_spacing_y,  # Adjust the spacing between legend items
        legend.margin = legend_margin #Use the custom legend margin
      ) + guides(
        color = guide_legend(
          override.aes = list(shape = 22, fill = colors_to_use, size = 6, color = "black")
        )
      ) + scale_shape_manual(values = rep(22, length(unique(seurat_object[[grouping_variable]]))))
      
      print(p)
    }
  })
  
  output$download_dimPlot <- downloadHandler(
    filename = function() {
      paste("dimension-plot", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      print(output$dimPlot())  # assuming output$dimPlot is a ggplot object
      dev.off()
    }
  )
  
  # Similar setup for PNG download, using ggsave from the ggplot2 package
  output$download_dimPlot_png <- downloadHandler(
    filename = function() {
      paste("dimension-plot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = output$dimPlot(), width = 10, height = 8, dpi = 300)
    }
  )
  
  # Gene Expression Plots
  output$genePlot <- renderPlot({
    req(input$geneSelection)  # Ensure a gene is selected
    
    gene <- input$geneSelection
    plotType <- input$plotType
    
    if(gene %in% rownames(seurat_object@assays$RNA@data)) {
      if(plotType == "Feature Plot") {
        # Feature Plot
        featurePlot <- FeaturePlot(seurat_object, pt.size = 0.75, min.cutoff = "q5", max.cutoff = "q95", features = gene)
        print(featurePlot)
      } else if(plotType == "Violin Plot") {
        # Violin Plot
        VlnPlot(seurat_object, features = gene, pt.size = 0, group.by = "Celltype",cols = celltype_colors) + NoLegend()
      } else if(plotType == "Box and whisker plot") {
        # Box and Whisker Plot (Bar Graph)
        data <- FetchData(seurat_object, vars = c("ident", gene))
        ggplot(data, aes(x = ident, y = !!sym(gene))) +
          geom_boxplot() +
          theme_minimal() +
          labs(x = "Identity", y = "Expression", title = paste("Expression of", gene)) +
          coord_flip()  # Optional: Flip coordinates for horizontal layout
      }
    } else {
      plot.new()
      text(0.5, 0.5, paste("Gene", gene, "not found"), cex = 2)
    }
  })

  # Gene Expression Plots
  output$genePlot2 <- renderPlot({
    req(input$geneSelection2)  # Ensure a gene is selected
    
    gene <- input$geneSelection2
    plotType <- input$plotType2
    
    if(gene %in% rownames(seurat_object2@assays[["RNA"]])) {
      if(plotType == "Feature Plot") {
        # Feature Plot
        featurePlot <- FeaturePlot(seurat_object2, pt.size = 2, min.cutoff = "q5", max.cutoff = "q95", features = gene)
        print(featurePlot)
      } else if(plotType == "Violin Plot") {
        # Violin Plot
        VlnPlot(seurat_object2, features = gene, pt.size = 0, group.by = "Celltype_Stage",cols = ordered_colors3) + NoLegend()
      } else if(plotType == "Box and whisker plot") {
        # Box and Whisker Plot (Bar Graph)
        data <- FetchData(seurat_object2, vars = c("ident", gene))
        ggplot(data, aes(x = ident, y = !!sym(gene))) +
          geom_boxplot() +
          theme_minimal() +
          labs(x = "Identity", y = "Expression", title = paste("Expression of", gene)) +
          coord_flip()  # Optional: Flip coordinates for horizontal layout
      }
    } else {
      plot.new()
      text(0.5, 0.5, paste("Gene", gene, "not found"), cex = 2)
    }
  })
  
  output$cellTypeDimPlot <- renderPlot({
    # Use the 'Celltype' grouping for this static DimPlot
    legend_margin <- margin(t = 0, r = 250, b = 0, l = 10, unit = "pt")
    p <- DimPlot(
      object = seurat_object,
      reduction = "umap",
      group.by = "Celltype",
      cols = ordered_colors,  # Assuming 'ordered_colors' maps correctly to 'Celltype'
      pt.size = 0.75
    ) + theme(
      legend.position = "right",
      legend.margin = legend_margin
    ) + guides(
      color = guide_legend(
        override.aes = list(shape = 22, fill = ordered_colors, size = 6, color = "black")
      )
    ) + scale_shape_manual(values = rep(22, length(unique(seurat_object$Celltype))))
    
    print(p)
  })

  output$cellTypeDimPlot2 <- renderPlot({
    # Use the 'Celltype_Stage' grouping for this static DimPlot
    legend_margin <- margin(t = 0, r = 250, b = 0, l = 10, unit = "pt")
    p <- DimPlot(
      object = seurat_object2,
      reduction = "umap",
      group.by = "Celltype_Stage",
      cols = ordered_colors3,  
      pt.size = 2
    ) + theme(
      legend.position = "right",
      legend.margin = legend_margin
    ) + guides(
      color = guide_legend(
        override.aes = list(shape = 22, fill = ordered_colors3, size = 6, color = "black")
      )
    ) + scale_shape_manual(values = rep(22, length(unique(seurat_object2$Celltype_Stage))))
    
    print(p)
  })
  
}

# Run the app
shinyApp(ui = ui, server = server)





