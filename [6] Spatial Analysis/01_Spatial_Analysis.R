library(SingleCellExperiment)
library(imcRtools)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggpubr)
library(viridis)
library(ggpointdensity)
library(singleCellTK)
library(pheatmap)
library(BiocParallel)
library(cytomapper)
library(scales)


#--------------------------------Loading Images------------------------------------#
images <- loadImages("./path/to/data",)
masks <- loadImages("./path/to/data",)


markers_columns <- c("Nuclear","CD31","CK5","SMA","Ki67","CK8","ANAX1","ECAD","PCNT","CCASP3")
class_colors <- c(Luminal = "#51abcb", BnL = "#f8766d",Basal = "red", `Blood Vessel` = "#a3a500", LnB = "#00b0f6",Undefined = "grey")
cells <- read.csv("./path/to/data")

counts <- cells[,markers_columns]
cell_meta <- cells[,!(names(cells) %in% markers_columns)]
sce <- SingleCellExperiment(assay = list(counts = t(counts)))
assay(sce,"exprs") <- as.matrix(t(counts))
colData(sce) <- as(cell_meta, "DataFrame")
spe <- buildSpatialGraph(sce, img_id = "imageNb", type = "expansion", threshold = 60)

spe <- readRDS("./path/to/data")

plotSpatial(spe[,spe$imageNb == "Mx_BEGIN TMA NORMAL_A07.ome_896_4032_1408_4544.tiff"],
            node_color_by = "class", 
            draw_edges = TRUE,
            colPairName = "expansion_interaction_graph", 
            img_id = "imageNb", 
            node_size_fix = 10)

findBorderCells(spe,
                img_id = "imageNb",
                border_dist = 10)

spe <- aggregateNeighbors(spe, colPairName = "expansion_interaction_graph", 
                          aggregate_by = "metadata", count_by = "class")

spe <- detectSpatialContext(spe, entry = "aggregatedNeighbors",
                            threshold = 0.9)

plotSpatialContext(spe, group_by = "imageNb",
                   node_color_by = "n_group",
                   node_size_by = "n_cells",
                   node_label_color_by = "n_group")+
  scale_color_identity()

spe <- patchDetection(spe, 
                      patch_cells = spe$class == "BnL",
                      img_id = "imageNb",
                      expand_by = 1,
                      min_patch_size = 10,
                      colPairName = "neighborhood",
                      BPPARAM = MulticoreParam())

cn_1 <- kmeans(spe$aggregatedNeighbors, centers = 8)
spe$cn_celltypes <- as.factor(cn_1$cluster)


plotSpatial(spe[,spe$imageNb == "Mx_BEGIN TMA EDGE_A10.ome_0_2240_512_2752.tiff"],
            node_color_by = "cn_celltypes",
            draw_edges = TRUE,
            colPairName = "expansion_interaction_graph", 
            img_id = "imageNb", 
            node_size_fix = 10)

plotSpatial(spe[,spe$imageNb == "Mx_BEGIN TMA NORMAL_A07.ome_448_2688_960_3200.tiff"],
            node_color_by = "cn_celltypes", 
            draw_edges = TRUE,
            colPairName = "expansion_interaction_graph", 
            img_id = "imageNb", 
            node_size_fix = 10)

plotSpatial(spe[,spe$imageNb == "Mx_BEGIN TMA NORMAL_A07.ome_448_2688_960_3200.tiff"],
            node_color_by = "class", 
            draw_edges = TRUE,
            colPairName = "expansion_interaction_graph", 
            img_id = "imageNb", 
            node_size_fix = 10)

# Heatmap of cell types in each neighborhood for PDG
for_plot <- colData(spe) %>% as_tibble() %>%
  group_by(cn_celltypes, class) %>%
  summarize(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  pivot_wider(id_cols = cn_celltypes, names_from = class, 
              values_from = freq, values_fill = 0) %>%
  ungroup() %>%
  select(-cn_celltypes)

pheatmap(for_plot,
         color=colorRampPalette(c("darkblue","white","darkred"))(100), 
         scale = "column")

saveRDS(spe,"./path/to/data")

#------------------------------------------------------------------------------------------------------------------#
cells_Normal <- cells[grepl("NORMAL", cells$imageNb, ignore.case = TRUE), ]
cells_Edge <- cells[grepl("EDGE", cells$imageNb, ignore.case = TRUE), ]
cells_Tumour <- cells[grepl("TUMOUR", cells$imageNb, ignore.case = TRUE), ]
cells_LuminalA <- cells[cells$diagnosis == "Luminal A",]
cells_LuminalB <- cells[cells$diagnosis == "Luminal B",]
cells_HER2 <- cells[cells$diagnosis == "HER2 enriched",]
cells_Triple <- cells[cells$diagnosis == "Triple negative",]
cells_DCIS <- cells[cells$Histology_description == "DCIS",]
cells_Invasive <- cells[cells$Histology_description == "Invasive ductal",]
cells_Lobular <- cells[cells$Histology_description == "Lobular",]
cells_Papillary <- cells[cells$Histology_description == "Papillary",]
cells_Mixed <- cells[cells$Histology_description == "Mixed",]
cells_Grade0 <- cells[cells$GRADE == 0,]
cells_Grade1 <- cells[cells$GRADE == 1,]
cells_Grade2 <- cells[cells$GRADE == 2,]
cells_Grade3 <- cells[cells$GRADE == 3,]

# Create SCE objects for each treatment condition
counts1 <- cells_Normal[,markers_columns]
cell_meta1 <- cells_Normal[,!(names(cells_Normal) %in% markers_columns)]
sce_Normal <- SingleCellExperiment(assays = list(counts1 = t(counts1)))
colData(sce_Normal) <- as(cell_meta1, "DataFrame")

counts2 <- cells_Edge[,markers_columns]
cell_meta2 <- cells_Edge[,!(names(cells_Edge) %in% markers_columns)]
sce_Edge <- SingleCellExperiment(assays = list(counts2 = t(counts2)))
colData(sce_Edge) <- as(cell_meta2, "DataFrame")

counts3 <- cells_Tumour[,markers_columns]
cell_meta3 <- cells_Tumour[,!(names(cells_Tumour) %in% markers_columns)]
sce_Tumour <- SingleCellExperiment(assays = list(counts3 = t(counts3)))
colData(sce_Tumour) <- as(cell_meta3, "DataFrame")

counts4 <- cells_LuminalA[,markers_columns]
cell_meta4 <- cells_LuminalA[,!(names(cells_LuminalA) %in% markers_columns)]
sce_LuminalA <- SingleCellExperiment(assays = list(counts4 = t(counts4)))
colData(sce_LuminalA) <- as(cell_meta4, "DataFrame")

counts5 <- cells_LuminalB[,markers_columns]
cell_meta5 <- cells_LuminalB[,!(names(cells_LuminalB) %in% markers_columns)]
sce_LuminalB <- SingleCellExperiment(assays = list(counts5 = t(counts5)))
colData(sce_LuminalB) <- as(cell_meta5, "DataFrame")

counts6 <- cells_HER2[,markers_columns]
cell_meta6 <- cells_HER2[,!(names(cells_HER2) %in% markers_columns)]
sce_HER2 <- SingleCellExperiment(assays = list(counts6 = t(counts6)))
colData(sce_HER2) <- as(cell_meta6, "DataFrame")

counts7 <- cells_Triple[,markers_columns]
cell_meta7 <- cells_Triple[,!(names(cells_Triple) %in% markers_columns)]
sce_Triple <- SingleCellExperiment(assays = list(counts7 = t(counts7)))
colData(sce_Triple) <- as(cell_meta7, "DataFrame")

counts8 <- cells_DCIS[,markers_columns]
cell_meta8 <- cells_DCIS[,!(names(cells_DCIS) %in% markers_columns)]
sce_DCIS <- SingleCellExperiment(assays = list(counts8 = t(counts8)))
colData(sce_DCIS) <- as(cell_meta8, "DataFrame")

counts9 <- cells_Invasive[,markers_columns]
cell_meta9 <- cells_Invasive[,!(names(cells_Invasive) %in% markers_columns)]
sce_Invasive <- SingleCellExperiment(assays = list(counts9 = t(counts9)))
colData(sce_Invasive) <- as(cell_meta9, "DataFrame")

counts10 <- cells_Lobular[,markers_columns]
cell_meta10 <- cells_Lobular[,!(names(cells_Lobular) %in% markers_columns)]
sce_Lobular <- SingleCellExperiment(assays = list(counts10 = t(counts10)))
colData(sce_Lobular) <- as(cell_meta10, "DataFrame")

counts11 <- cells_Papillary[,markers_columns]
cell_meta11 <- cells_Papillary[,!(names(cells_Papillary) %in% markers_columns)]
sce_Papillary <- SingleCellExperiment(assays = list(counts11 = t(counts11)))
colData(sce_Papillary) <- as(cell_meta11, "DataFrame")

counts12 <- cells_Mixed[,markers_columns]
cell_meta12 <- cells_Mixed[,!(names(cells_Mixed) %in% markers_columns)]
sce_Mixed <- SingleCellExperiment(assays = list(counts12 = t(counts12)))
colData(sce_Mixed) <- as(cell_meta12, "DataFrame")

counts13 <- cells_Grade0[,markers_columns]
cell_meta13 <- cells_Grade0[,!(names(cells_Grade0) %in% markers_columns)]
sce_Grade0 <- SingleCellExperiment(assays = list(counts13 = t(counts13)))
colData(sce_Grade0) <- as(cell_meta13, "DataFrame")

counts14 <- cells_Grade1[,markers_columns]
cell_meta14 <- cells_Grade1[,!(names(cells_Grade1) %in% markers_columns)]
sce_Grade1 <- SingleCellExperiment(assays = list(counts14 = t(counts14)))
colData(sce_Grade1) <- as(cell_meta14, "DataFrame")

counts15 <- cells_Grade2[,markers_columns]
cell_meta15 <- cells_Grade2[,!(names(cells_Grade2) %in% markers_columns)]
sce_Grade2 <- SingleCellExperiment(assays = list(counts15 = t(counts15)))
colData(sce_Grade2) <- as(cell_meta15, "DataFrame")

counts16 <- cells_Grade3[,markers_columns]
cell_meta16 <- cells_Grade3[,!(names(cells_Grade3) %in% markers_columns)]
sce_Grade3 <- SingleCellExperiment(assays = list(counts16 = t(counts16)))
colData(sce_Grade3) <- as(cell_meta16, "DataFrame")

# Create spatial interaction graph per image for each treatment condition using 'expansion' method set to 30 micron diamter. This step can take hours to a full day

spe_Normal <- buildSpatialGraph(sce_Normal, img_id = "imageNb", type = "expansion", threshold = 60)
spe_Edge <- buildSpatialGraph(sce_Edge, img_id = "imageNb", type = "expansion", threshold = 60)
spe_Tumour <- buildSpatialGraph(sce_Tumour, img_id = "imageNb", type = "expansion", threshold = 60)
spe_LuminalA <- buildSpatialGraph(sce_LuminalA, img_id = "imageNb", type = "expansion", threshold = 60)
spe_LuminalB <- buildSpatialGraph(sce_LuminalB, img_id = "imageNb", type = "expansion", threshold = 60)
spe_HER2 <- buildSpatialGraph(sce_HER2, img_id = "imageNb", type = "expansion", threshold = 60)
spe_Triple <- buildSpatialGraph(sce_Triple, img_id = "imageNb", type = "expansion", threshold = 60)
spe_DCIS <- buildSpatialGraph(sce_DCIS, img_id = "imageNb", type = "expansion", threshold = 60)
spe_Invasive <- buildSpatialGraph(sce_Invasive, img_id = "imageNb", type = "expansion", threshold = 60)
spe_Lobular <- buildSpatialGraph(sce_Lobular, img_id = "imageNb", type = "expansion", threshold = 60)
spe_Papillary <- buildSpatialGraph(sce_Papillary, img_id = "imageNb", type = "expansion", threshold = 60)
spe_Mixed <- buildSpatialGraph(sce_Mixed, img_id = "imageNb", type = "expansion", threshold = 60)
spe_Grade0 <- buildSpatialGraph(sce_Grade0, img_id = "imageNb", type = "expansion", threshold = 60)
spe_Grade1 <- buildSpatialGraph(sce_Grade1, img_id = "imageNb", type = "expansion", threshold = 60)
spe_Grade2 <- buildSpatialGraph(sce_Grade2, img_id = "imageNb", type = "expansion", threshold = 60)
spe_Grade3 <- buildSpatialGraph(sce_Grade3, img_id = "imageNb", type = "expansion", threshold = 60)

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# Function that aggregates and clusters cell types, followed by plotting onto a heatmap. 
generate_neighbours_heatmap <- function(obj,filename,colPairName = "expansion_interaction_graph", 
                                        aggregate_by = "metadata", count_by = "class", 
                                        centers = 6, color_palette = c("darkblue", "white", "darkred")) {
  
  # Aggregate neighbors
  obj <- aggregateNeighbors(obj, colPairName = colPairName, 
                            aggregate_by = aggregate_by, count_by = count_by)
  
  # Perform k-means clustering
  cn_1 <- kmeans(obj$aggregatedNeighbors, centers = centers)
  obj$cn_celltypes <- as.factor(cn_1$cluster)
  
  # Create heatmap data
  for_plot <- colData(obj) %>% 
    as_tibble() %>%
    group_by(cn_celltypes, class) %>%
    summarize(count = n(), .groups = 'drop') %>%
    mutate(freq = count / sum(count)) %>%
    pivot_wider(id_cols = cn_celltypes, names_from = class, 
                values_from = freq, values_fill = list(freq = 0)) %>%
    ungroup() %>%
    select(-cn_celltypes)
  
  # Generate the heatmap
  x <- pheatmap(for_plot,
           color = colorRampPalette(color_palette)(100), 
           scale = "column")
  
  save_pheatmap_pdf(x, filename)
  
  return(obj)
}

spe_Normal <- generate_neighbours_heatmap(obj = spe_Normal,filename = "./spe_Normal_freq.pdf")
spe_Edge <- generate_neighbours_heatmap(obj = spe_Edge,filename = "./spe_Edge_freq.pdf")
spe_Tumour <- generate_neighbours_heatmap(obj = spe_Tumour,filename = "./spe_Tumour_freq.pdf")
spe_LuminalA <- generate_neighbours_heatmap(obj = spe_LuminalA,filename = "./spe_LuminalA_freq.pdf")
spe_LuminalB <- generate_neighbours_heatmap(obj = spe_LuminalB,filename = "./spe_LuminalB_freq.pdf")
spe_HER2 <- generate_neighbours_heatmap(obj = spe_HER2,filename = "./spe_HER2_freq.pdf")
spe_Triple <- generate_neighbours_heatmap(obj = spe_Triple,filename = "./spe_Triple_freq.pdf")
spe_DCIS <- generate_neighbours_heatmap(obj = spe_DCIS,filename = "./spe_DCIS_freq.pdf")
spe_Invasive <- generate_neighbours_heatmap(obj = spe_Invasive,filename = "./spe_Invasive_freq.pdf")
spe_Lobular <- generate_neighbours_heatmap(obj = spe_Lobular,filename = "./spe_Lobular_freq.pdf")
spe_Papillary <- generate_neighbours_heatmap(obj = spe_Papillary,filename = "./spe_Papillary_freq.pdf")
spe_Mixed <- generate_neighbours_heatmap(obj = spe_Mixed,filename = "./spe_Mixed_freq.pdf")
spe_Grade0 <- generate_neighbours_heatmap(obj = spe_Grade0,filename = "./spe_Grade0_freq.pdf")
spe_Grade1 <- generate_neighbours_heatmap(obj = spe_Grade1,filename = "./spe_Grade1_freq.pdf")
spe_Grade2 <- generate_neighbours_heatmap(obj = spe_Grade2,filename = "./spe_Grade2_freq.pdf")
spe_Grade3 <- generate_neighbours_heatmap(obj = spe_Grade3,filename = "./spe_Grade3_freq.pdf")

saveRDS(spe_Normal, "./spe_Normal.rds")
saveRDS(spe_Edge, "./spe_Edge.rds")
saveRDS(spe_Tumour, "./spe_Tumour.rds")
saveRDS(spe_LuminalA, "./spe_LuminalA.rds")
saveRDS(spe_LuminalB, "./spe_LuminalB.rds")
saveRDS(spe_HER2, "./spe_HER2.rds")
saveRDS(spe_Triple, "./spe_Triple.rds")
saveRDS(spe_DCIS, "./spe_DCIS.rds")
saveRDS(spe_Invasive, "./spe_Invasive.rds")
saveRDS(spe_Lobular, "./spe_Lobular.rds")
saveRDS(spe_Papillary, "./spe_Papillary.rds")
saveRDS(spe_Mixed, "./spe_Mixed.rds")
saveRDS(spe_Grade0, "./spe_Grade0.rds")
saveRDS(spe_Grade1, "./spe_Grade1.rds")
saveRDS(spe_Grade2, "./spe_Grade2.rds")
saveRDS(spe_Grade3, "./spe_Grade3.rds")
