
# Install and load required packages
pkgs <- c("Seurat", "tidyverse", "gridExtra", "clustree", "sctransform", "patchwork",
          "RColorBrewer", "ggrepel", "viridis", "scCustomize", "ComplexHeatmap")

options(repos = c(CRAN = "https://cloud.r-project.org"))

missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing)) {
  cran <- intersect(missing, rownames(available.packages()))
  bioc <- setdiff(missing, cran)
  if (length(cran)) install.packages(cran, quiet = TRUE)
  if (length(bioc)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", quiet = TRUE)
    BiocManager::install(bioc, ask = FALSE, update = FALSE, quiet = TRUE)
  }
}

invisible(lapply(pkgs, library, character.only = TRUE))

# Set up environment
dir.create("~/human_otic/")
setwd("~/human_otic/")
getwd()
dir.create("plots")
dir.create("rds_files")
plot_path = "./plots/"
rds_path = "./rds_files/"

# Load dataset from 

Day0_36_all <- readRDS(paste0(rds_path,"GSE188936_Day0-Day36.rds"))

# Convert active.ident to cell.type
Day0_36_all@meta.data$celltype <- Day0_36_all@active.ident

# Calculate % Ribosomal genes
Day0_36_all[["percent.Ribosomal"]] <- PercentageFeatureSet(
  Day0_36_all,
  pattern = "^RP[LS]"
)

# Print Day (sample) and Cell Type Dimension Plots

d1 <- DimPlot(Day0_36_all, reduction = "umap", label = TRUE) + 
  theme(aspect.ratio = 1) +
  labs(title = "Human Otic cells\nDay0 to 36\noriginal assignments")

Day0_36_all_dimplot_celltype

ggsave(paste0(plot_path, "Day0_36_all_dimplot_celltype.png"), 
       d1, 
       units = "cm", height = 24, width = 24)

d2 <- DimPlot(Day0_36_all, 
              group.by = "sample",
              reduction = "umap", label = FALSE) + 
  theme(aspect.ratio = 1) +
  labs(title = "Human Otic cells\nDay0 to 36")


ggsave(paste0(plot_path, "Day0_36_all_dimplot_day.png"), 
       d2, 
       units = "cm", height = 24, width = 24)


#############################
# Subset Otic cells from Days 18 to 36

stages <- c("day 18", "day 24", "day 30", "day 36")

Otic_Day18_36 <- subset(Day0_36_all,
                   subset = sample %in% stages & celltype == "Otic")

#########################
# Recluster using SCTransform

Otic_Day18_36 <- SCTransform(Otic_Day18_36,
                             ncells = ncol(Otic_Day18_36),
                             verbose = F
) 

Otic_Day18_36 <- RunPCA(Otic_Day18_36, verbose = F)

Otic_Day18_36 <- FindNeighbors(Otic_Day18_36, dims = 1:30, verbose = F, reduction = "pca")
Otic_Day18_36 <- FindClusters(Otic_Day18_36, verbose = F, cluster.name = "clst.Otic_Day18_36")
Otic_Day18_36 <- RunUMAP(Otic_Day18_36, dims = 1:30, verbose = F, reduction = "pca", reduction.name = "umap.otic")

######################
# Testing resolution
# ==============================================================================
# CLUSTERING RESOLUTION SWEEP
# ==============================================================================
res_range <- seq(0.4, 2.0, by = 0.2)

for (res in res_range) {
  col_name <- paste0("otic_cluster_res", res)
  
  Otic_Day18_36 <- FindClusters(
    object       = Otic_Day18_36,
    resolution   = res,
    cluster.name = col_name,
    verbose      = FALSE
  )
  
  # DimPlot for each resolution
  plot_name <- paste0("Otic_Day18_36_dimplot_clusterres", res)
  p <- DimPlot(Otic_Day18_36,
               group.by  = col_name,
               reduction = "umap.otic",
               label     = TRUE,
               repel     = TRUE) +
    theme(aspect.ratio = 2) +
    labs(title = paste0("Resolution: ", res))
  
  assign(plot_name, p)
  
  ggsave(paste0(plot_path, plot_name, ".png"),
         p,
         units  = "cm",
         height = 24,
         width  = 14)
}

# ==============================================================================
# CLUSTREE
# ==============================================================================
tree_plot <- clustree(Otic_Day18_36, prefix = "otic_cluster_res")
tree_plot

ggsave(
  filename = file.path(plot_path, "Otic_Day18_36_clustree.png"),
  plot     = tree_plot,
  width    = 24,
  height   = 24,
  units    = "cm",
  dpi      = 300
)


###############################
# From clustree optimal resolution = 1.2 (otic_cluster_res1.2)
# 
# Print Day (sample) and Cluster Dimension Plots

# Cluster
cluster_cols <- c(
  "0"  = "#A1005E",
  "1"  = "#524C50",
  "2"  = "#F400ED",
  "3"  = "#2977FF",
  "4"  = "#24D0FF",
  "5"  = "#00F916",
  "6"  = "#8F00FF",
  "7"  = "#ED2025",
  "8"  = "#FA9E11",
  "9"  = "#D993FF",
  "10" = "#FB9D97",
  "11" = "#27FFCB",
  "12" = "#DFDDE2",
  "13" = "#8BAC17"
)

d3 <- DimPlot(Otic_Day18_36, 
      group.by = "otic_cluster_res1.2",
      reduction = "umap.otic", 
      label = TRUE,
      label.size = 8,
      cols = cluster_cols,
      pt.size = 2) + 
  theme(aspect.ratio = 2) +
  labs(title = NULL,
       x = NULL,
       y = NULL)

ggsave(paste0(plot_path, 'Otic_Day18_36_cluster_res1.2_recoloured.png'),
       d3,
       height = 24, width = 24, units = 'cm', dpi = 400)

# Day
d4 <- DimPlot(Otic_Day18_36, 
              group.by = "sample",
              reduction = "umap.otic", 
              label = F,
              pt.size = 2) + 
  theme(aspect.ratio = 2) +
  labs(title = NULL,
       x = NULL,
       y = NULL)

ggsave(paste0(plot_path, 'Otic_Day18_36_cluster_Day.png'),
       d4,
       height = 24, width = 24, units = 'cm', dpi = 400)

#####################
##################
# QC plots for res 1.2

v1 <- VlnPlot_scCustom(Otic_Day18_36, features = "nCount_RNA", group.by = "otic_cluster_res1.2")
v2 <- VlnPlot_scCustom(Otic_Day18_36, features = "nFeature_RNA", group.by = "otic_cluster_res1.2")
v3 <- VlnPlot_scCustom(Otic_Day18_36, features = "percent.Ribosomal", group.by = "otic_cluster_res1.2")
v4 <- VlnPlot_scCustom(Otic_Day18_36, features = "log10GenesPerUMI", group.by = "otic_cluster_res1.2")

ggsave(paste0(plot_path, "Otic_Day18_36_vlnplot_nCount_RNA.png"), v1, 
       units = "cm", height = 10, width = 14, dpi = 600)

ggsave(paste0(plot_path, "Otic_Day18_36_vlnplot_nFeature_RNA.png"), v2, 
       units = "cm", height = 10, width = 14, dpi = 600)

ggsave(paste0(plot_path, "Otic_Day18_36_vlnplot_percent.Ribosomal.png"), v3, 
       units = "cm", height = 10, width = 14, dpi = 600)

ggsave(paste0(plot_path, "Otic_Day18_36_vlnplot_log10GenesPerUMI.png"), v4, 
       units = "cm", height = 10, width = 14, dpi = 600)

######################################
# FindallMarkers to look at cluster enriched genes.

Idents(Otic_Day18_36) <- "otic_cluster_res1.2"

all_markers <- FindAllMarkers(object = Otic_Day18_36)%>%
  Add_Pct_Diff()

write_csv(all_markers, "Otic_Day18_36_allmarkers_clst1.2.csv")

all_markers_06 <- all_markers %>%
  filter(abs(pct_diff) > 0.6)

#############################################
# Feature plots for selected genes

# Function to construct feature plots for selected genes.

fp_otic <- function(sel.gene, stage, reduct){
  
  result <- tryCatch({
    
    fp1 <- FeaturePlot_scCustom(
      stage,
      features = sel.gene,
      reduction = reduct
    ) + theme(aspect.ratio = 1)
    
  }, error = function(e) {
    cat("Error occurred:", conditionMessage(e), "\n")
    return(NULL)
  }, warning = function(w) {
    cat("Warning occurred:", conditionMessage(w), "\n")
    return(NULL)
  })
  
  if (inherits(result, "try-error")) {
    return(NULL)
  }
  
  
  fp1 <- FeaturePlot_scCustom(
    stage,
    features = sel.gene,
    reduction = reduct,
    pt.size = 1,
    order = TRUE
  ) +
    theme(aspect.ratio = 2)
  
  # --------------------------------------------------------
  
  dir_name <- deparse(substitute(stage))
  dir.create(paste0("./feature_plots/", dir_name), showWarnings = FALSE)
  
  plot_path <- paste0("./feature_plots/", dir_name, "/")
  
  fp.filename <- paste0(sel.gene, "_", dir_name, "_featureplot_sc.png")
  ggsave(
    paste0(plot_path, fp.filename),
    fp1,
    height = 12,
    width = 10,
    units = "cm",
    dpi = 400
  )
  
  cat("Plot saved as", paste0(plot_path, fp.filename), "\n")
  
  grid.arrange(fp1)
}

# Feature plots & Dotplots

Marker_genes <- c("EYA1", "PAX2", "IRX2", "ELAVL2", "ELAVL3", "ELAVL4", "POU4F1", 
                  "NEUROD1", "PRRX1", "PDGFRA", "COL1A1", "COL1A2", "EMP3", "ZEB1", 
                  "EMX2", "DACH1", "LMX1A", "WNT2B", "NOX3", "PAX8", "SPP1", "DLX5", 
                  "DLX6", "MSX1", "OTOG", "FGF20", "OTOL1", "USH1C", "ISL1", "FGF10", 
                  "JAG1", "BIRC5", "TOP2A", "MKI67", "UBE2C", "CDC20", "PLK1", 
                  "CCNB1", "MCM5", "MCM2", "MCM6", "CDC6", "TK1", "CDC45", "MCM4", 
                  "ALDH1A1", "EBF1", "OC90", "BMP4", "ID2", "ID3")

lapply(Marker_genes, fp_otic, Otic_Day18_36, "umap.otic")


dot1 <- Clustered_DotPlot(seurat_object = Otic_Day18_36, 
                        features = Marker_genes,
                        group.by = "cluster.otic_res1.2",
                        k = 9)

png(filename = paste0(plot_path,"Otic_Day18_36_clustdotplot.png"), 
    height = 18,
    width = 16,
    units = "cm",
    res = 800
)
draw(dot1[[2]])
dev.off()

######################
# Blended Dimplots

# 1. Define the gene pairs using a tibble (tidy data frame)
gene_pairs <- tribble(
  ~gene1,  ~gene2,
  "WNT2B", "PAX2",
  "WNT2B", "DLX5",  "WNT2B", "HMX3",  "WNT2B", "PAX8", 
  "WNT2B", "JAG1",  "WNT2B", "FGF3",  "WNT2B", "FGF10",
  "DLX5",  "HMX3",  "DLX5",  "PAX8",  "DLX5",  "PAX2", 
  "DLX5",  "JAG1",  "DLX5",  "FGF3",  "DLX5",  "FGF10",
  "HMX3",  "PAX8",  "HMX3",  "PAX2",  "HMX3",  "JAG1", 
  "HMX3",  "FGF3",  "HMX3",  "FGF10", "PAX8",  "PAX2",
  "PAX8",  "JAG1",  "PAX8",  "FGF3",  "PAX8",  "FGF10",
  "PAX2",  "JAG1",  "PAX2",  "FGF3",  "PAX2",  "FGF10",
  "JAG1",  "FGF3",  "JAG1",  "FGF10", "FGF3",  "FGF10"
)

# 2. Use walk2 to iterate and save
walk2(gene_pairs$gene1, gene_pairs$gene2, function(g1, g2) {
  
  # Construct plot name dynamically
  plot_name <- str_glue("Otic_Day36_density_{g1}_{g2}")
  
  # Generate plot
  p <- Plot_Density_Joint_Only(
    seurat_object = Otic_Day36, 
    features = c(g1, g2), 
    reduction = "umap.otic"
  ) + 
    theme(aspect.ratio = 2) +
    theme(axis.title = element_blank())
  
  
  # Save plot using ggsave
  ggsave(
    filename = str_glue("{plot_path}{plot_name}.png"),
    plot = p,
    height = 12,
    width = 10,
    units = "cm",
    dpi = 400
  )
  
  # Optional: print progress to console
  cli::cli_alert_success("Saved: {plot_name}")
})


##################################################################
# Day 36 Only analysis
# Subset just Day36 but keep same umap.
#################################################################

Otic_Day36 <- subset(Otic_Day18_36,
                        subset = sample == "Day 36")


Otic_Day36 <- SCTransform(Otic_Day36,
                             ncells = ncol(Otic_Day36),
                             verbose = F) 

#########################################
# Print Cluster Dimension Plot

# Cluster
cluster_cols <- c(
  "0"  = "#A1005E",
  "1"  = "#524C50",
  "2"  = "#F400ED",
  "3"  = "#2977FF",
  "4"  = "#24D0FF",
  "5"  = "#00F916",
  "6"  = "#8F00FF",
  "7"  = "#ED2025",
  "8"  = "#FA9E11",
  "9"  = "#D993FF",
  "10" = "#FB9D97",
  "11" = "#27FFCB",
  "12" = "#DFDDE2",
  "13" = "#8BAC17"
)

d5 <- DimPlot(Otic_Day36, 
              group.by = "otic_cluster_res1.2",
              reduction = "umap.otic", 
              label = TRUE,
              label.size = 8,
              cols = cluster_cols,
              pt.size = 2) + 
  theme(aspect.ratio = 2) +
  labs(title = NULL,
       x = NULL,
       y = NULL)

ggsave(paste0(plot_path, 'Otic_Day36_cluster_res1.2_recoloured.png'),
       d5,
       height = 24, width = 24, units = 'cm', dpi = 400)


#####################
##################
# QC plots for res 1.2

v5 <- VlnPlot_scCustom(Otic_Day36, features = "nCount_RNA", group.by = "otic_cluster_res1.2")
v6 <- VlnPlot_scCustom(Otic_Day36, features = "nFeature_RNA", group.by = "otic_cluster_res1.2")
v7 <- VlnPlot_scCustom(Otic_Day36, features = "percent.Ribosomal", group.by = "otic_cluster_res1.2")
v8 <- VlnPlot_scCustom(Otic_Day36, features = "log10GenesPerUMI", group.by = "otic_cluster_res1.2")

ggsave(paste0(plot_path, "Otic_Day36_vlnplot_nCount_RNA.png"), v5, 
       units = "cm", height = 10, width = 14, dpi = 600)

ggsave(paste0(plot_path, "Otic_Day36_vlnplot_nFeature_RNA.png"), v6, 
       units = "cm", height = 10, width = 14, dpi = 600)

ggsave(paste0(plot_path, "Otic_Day36_vlnplot_percent.Ribosomal.png"), v7, 
       units = "cm", height = 10, width = 14, dpi = 600)

ggsave(paste0(plot_path, "Otic_Day36_vlnplot_log10GenesPerUMI.png"), v8, 
       units = "cm", height = 10, width = 14, dpi = 600)

######################################
# FindallMarkers to look at cluster enriched genes.

Idents(Otic_Day36) <- "otic_cluster_res1.2"

all_markers <- FindAllMarkers(object = Otic_Day36)%>%
  Add_Pct_Diff()

write_csv(all_markers, "Otic_Day36_allmarkers_clst1.2.csv")

all_markers_06 <- all_markers %>%
  filter(abs(pct_diff) > 0.6)

#############################################
# Feature plots for selected genes

# Function to construct feature plots for selected genes.

fp_otic <- function(sel.gene, stage, reduct){
  
  result <- tryCatch({
    
    fp1 <- FeaturePlot_scCustom(
      stage,
      features = sel.gene,
      reduction = reduct
    ) + theme(aspect.ratio = 1)
    
  }, error = function(e) {
    cat("Error occurred:", conditionMessage(e), "\n")
    return(NULL)
  }, warning = function(w) {
    cat("Warning occurred:", conditionMessage(w), "\n")
    return(NULL)
  })
  
  if (inherits(result, "try-error")) {
    return(NULL)
  }
  
  
  fp1 <- FeaturePlot_scCustom(
    stage,
    features = sel.gene,
    reduction = reduct,
    pt.size = 1,
    order = TRUE
  ) +
    theme(aspect.ratio = 2)
  
  # --------------------------------------------------------
  
  dir_name <- deparse(substitute(stage))
  dir.create(paste0("./feature_plots/", dir_name), showWarnings = FALSE)
  
  plot_path <- paste0("./feature_plots/", dir_name, "/")
  
  fp.filename <- paste0(sel.gene, "_", dir_name, "_featureplot_sc.png")
  ggsave(
    paste0(plot_path, fp.filename),
    fp1,
    height = 12,
    width = 10,
    units = "cm",
    dpi = 400
  )
  
  cat("Plot saved as", paste0(plot_path, fp.filename), "\n")
  
  grid.arrange(fp1)
}

# Feature plots & Dotplots

Marker_genes <- c("EYA1", "PAX2", "IRX2", "ELAVL2", "ELAVL3", "ELAVL4", "POU4F1", 
                  "NEUROD1", "PRRX1", "PDGFRA", "COL1A1", "COL1A2", "EMP3", "ZEB1", 
                  "EMX2", "DACH1", "LMX1A", "WNT2B", "NOX3", "PAX8", "SPP1", "DLX5", 
                  "DLX6", "MSX1", "OTOG", "FGF20", "OTOL1", "USH1C", "ISL1", "FGF10", 
                  "JAG1", "BIRC5", "TOP2A", "MKI67", "UBE2C", "CDC20", "PLK1", 
                  "CCNB1", "MCM5", "MCM2", "MCM6", "CDC6", "TK1", "CDC45", "MCM4", 
                  "ALDH1A1", "EBF1", "OC90", "BMP4", "ID2", "ID3")

lapply(Marker_genes, fp_otic, Otic_Day36, "umap.otic")


dot2 <- Clustered_DotPlot(seurat_object = Otic_Day36, 
                        features = Marker_genes,
                        group.by = "cluster.otic_res1.2",
                        k = 9)

png(filename = paste0(plot_path,"Otic_Day36_clustdotplot.png"), 
    height = 18,
    width = 16,
    units = "cm",
    res = 800
)
draw(dot2[[2]])
dev.off()

######################
# Blended Dimplots

# 1. Define the gene pairs using a tibble (tidy data frame)
gene_pairs <- tribble(
  ~gene1,  ~gene2,
  "WNT2B", "PAX2",
  "WNT2B", "DLX5",  "WNT2B", "HMX3",  "WNT2B", "PAX8", 
  "WNT2B", "JAG1",  "WNT2B", "FGF3",  "WNT2B", "FGF10",
  "DLX5",  "HMX3",  "DLX5",  "PAX8",  "DLX5",  "PAX2", 
  "DLX5",  "JAG1",  "DLX5",  "FGF3",  "DLX5",  "FGF10",
  "HMX3",  "PAX8",  "HMX3",  "PAX2",  "HMX3",  "JAG1", 
  "HMX3",  "FGF3",  "HMX3",  "FGF10", "PAX8",  "PAX2",
  "PAX8",  "JAG1",  "PAX8",  "FGF3",  "PAX8",  "FGF10",
  "PAX2",  "JAG1",  "PAX2",  "FGF3",  "PAX2",  "FGF10",
  "JAG1",  "FGF3",  "JAG1",  "FGF10", "FGF3",  "FGF10"
)

# 2. Use walk2 to iterate and save
walk2(gene_pairs$gene1, gene_pairs$gene2, function(g1, g2) {
  
  # Construct plot name dynamically
  plot_name <- str_glue("Otic_Day36_density_{g1}_{g2}")
  
  # Generate plot
  p <- Plot_Density_Joint_Only(
    seurat_object = Otic_Day36, 
    features = c(g1, g2), 
    reduction = "umap.otic"
  ) + 
    theme(aspect.ratio = 2) +
    theme(axis.title = element_blank())
  
  
  # Save plot using ggsave
  ggsave(
    filename = str_glue("{plot_path}{plot_name}.png"),
    plot = p,
    height = 12,
    width = 10,
    units = "cm",
    dpi = 400
  )
  
  # Print progress to console
  cli::cli_alert_success("Saved: {plot_name}")
})


