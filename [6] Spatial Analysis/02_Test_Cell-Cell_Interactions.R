
spe_Normal <- readRDS("E:/Final Data/spe_Normal.rds")
spe_Edge <- readRDS("E:/Final Data/spe_Edge.rds")
spe_Tumour <- readRDS("E:/Final Data/spe_Tumour.rds")
spe_LuminalA <- readRDS("E:/Final Data/spe_LuminalA.rds")
spe_LuminalB <- readRDS("E:/Final Data/spe_LuminalB.rds")
spe_HER2 <- readRDS("E:/Final Data/spe_HER2.rds")
spe_Triple <- readRDS("E:/Final Data/spe_Triple.rds")
spe_DCIS <- readRDS("E:/Final Data/spe_DCIS.rds")
spe_Invasive <- readRDS("E:/Final Data/spe_Invasive.rds")
spe_Lobular <- readRDS("E:/Final Data/spe_Lobular.rds")
spe_Papillary <- readRDS("E:/Final Data/spe_Papillary.rds")
spe_Mixed <- readRDS("E:/Final Data/spe_Mixed.rds")
spe_Grade0 <- readRDS("E:/Final Data/spe_Grade0.rds")
spe_Grade1 <- readRDS("E:/Final Data/spe_Grade1.rds")
spe_Grade2 <- readRDS("E:/Final Data/spe_Grade2.rds")
spe_Grade3 <- readRDS("E:/Final Data/spe_Grade3.rds")

# 'Classic' method used to test interactions between cell types. Can take multiple days to complete each of the following steps
# Heatmaps of cellular colocalization calculated as the sum of two one-tailed permutation test p values (sum_sigval) for (i) untreated samples and (j) 7 days post-IR samples.
out_Normal <- testInteractions(spe_Normal, 
                               group_by = "imageNb",
                               label = "class", 
                               iter = 1000,
                               colPairName = "expansion_interaction_graph",
                               BPPARAM = SerialParam(RNGseed = 123))
head(out_Normal)

out_Edge <- testInteractions(spe_Edge, 
                             group_by = "imageNb",
                             label = "class", 
                             colPairName = "expansion_interaction_graph",
                             iter = 1000,
                             BPPARAM = SerialParam(RNGseed = 123))
head(out_Edge)

out_Tumour <- testInteractions(spe_Tumour, 
                               group_by = "imageNb",
                               label = "class", 
                               colPairName = "expansion_interaction_graph",
                               iter = 1000,
                               BPPARAM = SerialParam(RNGseed = 123))
head(out_Tumour)

out_LuminalA <- testInteractions(spe_LuminalA, 
                               group_by = "imageNb",
                               label = "class", 
                               colPairName = "expansion_interaction_graph",
                               iter = 1000,
                               BPPARAM = SerialParam(RNGseed = 123))
head(out_LuminalA)

out_LuminalB <- testInteractions(spe_LuminalB, 
                                 group_by = "imageNb",
                                 label = "class", 
                                 colPairName = "expansion_interaction_graph",
                                 iter = 1000,
                                 BPPARAM = SerialParam(RNGseed = 123))
head(out_LuminalB)

out_HER2 <- testInteractions(spe_HER2, 
                                 group_by = "imageNb",
                                 label = "class", 
                                 colPairName = "expansion_interaction_graph",
                                 iter = 1000,
                                 BPPARAM = SerialParam(RNGseed = 123))
head(out_HER2)

out_Triple <- testInteractions(spe_Triple, 
                                 group_by = "imageNb",
                                 label = "class", 
                                 colPairName = "expansion_interaction_graph",
                                 iter = 1000,
                                 BPPARAM = SerialParam(RNGseed = 123))
head(out_Triple)

out_DCIS <- testInteractions(spe_DCIS, 
                                 group_by = "imageNb",
                                 label = "class", 
                                 colPairName = "expansion_interaction_graph",
                                 iter = 1000,
                                 BPPARAM = SerialParam(RNGseed = 123))
head(out_DCIS)

out_Invasive <- testInteractions(spe_Invasive, 
                                group_by = "imageNb",
                                label = "class", 
                                colPairName = "expansion_interaction_graph",
                                iter = 1000,
                                BPPARAM = SerialParam(RNGseed = 123))
head(out_Invasive)

out_Lobular <- testInteractions(spe_Lobular, 
                                 group_by = "imageNb",
                                 label = "class", 
                                 colPairName = "expansion_interaction_graph",
                                 iter = 1000,
                                 BPPARAM = SerialParam(RNGseed = 123))
head(out_Lobular)

out_Papillary <- testInteractions(spe_Papillary, 
                                 group_by = "imageNb",
                                 label = "class", 
                                 colPairName = "expansion_interaction_graph",
                                 iter = 1000,
                                 BPPARAM = SerialParam(RNGseed = 123))
head(out_Papillary)

out_Mixed <- testInteractions(spe_Mixed, 
                                  group_by = "imageNb",
                                  label = "class", 
                                  colPairName = "expansion_interaction_graph",
                                  iter = 1000,
                                  BPPARAM = SerialParam(RNGseed = 123))
head(out_Papillary)

out_Grade0 <- testInteractions(spe_Grade0, 
                              group_by = "imageNb",
                              label = "class", 
                              colPairName = "expansion_interaction_graph",
                              iter = 1000,
                              BPPARAM = SerialParam(RNGseed = 123))
head(out_Grade0)

out_Grade1 <- testInteractions(spe_Grade1, 
                               group_by = "imageNb",
                               label = "class", 
                               colPairName = "expansion_interaction_graph",
                               iter = 1000,
                               BPPARAM = SerialParam(RNGseed = 123))
head(out_Grade1)

out_Grade2 <- testInteractions(spe_Grade2, 
                               group_by = "imageNb",
                               label = "class", 
                               colPairName = "expansion_interaction_graph",
                               iter = 1000,
                               BPPARAM = SerialParam(RNGseed = 123))
head(out_Grade2)

out_Grade3 <- testInteractions(spe_Grade3, 
                               group_by = "imageNb",
                               label = "class", 
                               colPairName = "expansion_interaction_graph",
                               iter = 1000,
                               BPPARAM = SerialParam(RNGseed = 123))
head(out_Grade3)

write.csv(as.data.frame(out_Normal),"E:/Final Data/Normal_Interactions.csv",row.names = FALSE)
write.csv(as.data.frame(out_Edge),"E:/Final Data/Edge_Interactions.csv",row.names = FALSE)
write.csv(as.data.frame(out_Tumour),"E:/Final Data/Tumour_Interactions.csv",row.names = FALSE)
write.csv(as.data.frame(out_LuminalA),"E:/Final Data/LuminalA_Interactions.csv",row.names = FALSE)
write.csv(as.data.frame(out_LuminalB),"E:/Final Data/LuminalB_Interactions.csv",row.names = FALSE)
write.csv(as.data.frame(out_HER2),"E:/Final Data/HER2_Interactions.csv",row.names = FALSE)
write.csv(as.data.frame(out_Triple),"E:/Final Data/Triple_Interactions.csv",row.names = FALSE)
write.csv(as.data.frame(out_DCIS),"E:/Final Data/DCIS_Interactions.csv",row.names = FALSE)
write.csv(as.data.frame(out_Invasive),"E:/Final Data/Invasive_Interactions.csv",row.names = FALSE)
write.csv(as.data.frame(out_Lobular),"E:/Final Data/Lobular_Interactions.csv")
write.csv(as.data.frame(out_Papillary),"E:/Final Data/Papillary_Interactions.csv",row.names = FALSE)
write.csv(as.data.frame(out_Mixed),"E:/Final Data/Mixed_Interactions.csv",row.names = FALSE)
write.csv(as.data.frame(out_Grade0),"E:/Final Data/Grade0_Interactions.csv",row.names = FALSE)
write.csv(as.data.frame(out_Grade1),"E:/Final Data/Grade1_Interactions.csv",row.names = FALSE)
write.csv(as.data.frame(out_Grade2),"E:/Final Data/Grade2_Interactions.csv",row.names = FALSE)
write.csv(as.data.frame(out_Grade3),"E:/Final Data/Grade3_Interactions.csv",row.names = FALSE)


out_Normal <- read.csv("E:/Final Data/Normal_Interactions.csv", header = TRUE)
out_Edge <- read.csv("E:/Final Data/Edge_Interactions.csv", header = TRUE)
out_Tumour <- read.csv("E:/Final Data/Tumour_Interactions.csv", header = TRUE)
out_LuminalA <- read.csv("E:/Final Data/LuminalA_Interactions.csv", header = TRUE)
out_LuminalB <- read.csv("E:/Final Data/LuminalB_Interactions.csv", header = TRUE)
out_HER2 <- read.csv("E:/Final Data/HER2_Interactions.csv", header = TRUE)
out_Triple <- read.csv("E:/Final Data/Triple_Interactions.csv", header = TRUE)
out_DCIS <- read.csv("E:/Final Data/DCIS_Interactions.csv", header = TRUE)
out_Invasive <- read.csv("E:/Final Data/Invasive_Interactions.csv", header = TRUE)
out_Lobular <- read.csv("E:/Final Data/Lobular_Interactions.csv", header = TRUE)
out_Papillary <- read.csv("E:/Final Data/Papillary_Interactions.csv", header = TRUE)
out_Mixed <- read.csv("E:/Final Data/Mixed_Interactions.csv", header = TRUE)
out_Grade0 <- read.csv("E:/Final Data/Grade0_Interactions.csv", header = TRUE)
out_Grade1 <- read.csv("E:/Final Data/Grade1_Interactions.csv", header = TRUE)
out_Grade2 <- read.csv("E:/Final Data/Grade2_Interactions.csv", header = TRUE)
out_Grade3 <- read.csv("E:/Final Data/Grade3_Interactions.csv", header = TRUE)


out_Normal %>% as_tibble() %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>% 
  mutate(log_sigval = log1p(abs(sum_sigval)) * sign(sum_sigval)) %>%  # Log transformation with sign preservation
  ggplot(aes(x = from_label, y = to_label, fill = log_sigval)) +
  geom_tile() +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"), midpoint = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

out_Edge %>% as_tibble() %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>% 
  mutate(log_sigval = log1p(abs(sum_sigval)) * sign(sum_sigval)) %>%  # Log transformation with sign preservation
  ggplot(aes(x = from_label, y = to_label, fill = log_sigval)) +
  geom_tile() +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"), midpoint = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

out_Tumour %>% as_tibble() %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>% 
  mutate(log_sigval = log1p(abs(sum_sigval)) * sign(sum_sigval)) %>%  # Log transformation with sign preservation
  ggplot(aes(x = from_label, y = to_label, fill = log_sigval)) +
  geom_tile() +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"), midpoint = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


heatmap_sigval <- function(df,fpath){
  plot <- df %>% as_tibble() %>%
    group_by(from_label, to_label) %>%
    summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>% 
    mutate(log_sigval = log1p(abs(sum_sigval)) * sign(sum_sigval)) %>%  # Log transformation with sign preservation
    ggplot(aes(x = from_label, y = to_label, fill = log_sigval)) +
    geom_tile() +
    scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"), midpoint = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  
    ggsave(fpath)
}

heatmap_sigval(out_Normal,"./normal_sum_sigval.pdf")
heatmap_sigval(out_Edge,"./edge_sum_sigval.pdf")
heatmap_sigval(out_Tumour,"./tumour_sum_sigval.pdf")
heatmap_sigval(out_LuminalA,"./luminalA_sum_sigval.pdf")
heatmap_sigval(out_LuminalB,"./luminalB_sum_sigval.pdf")
heatmap_sigval(out_HER2,"./HER2_sum_sigval.pdf")
heatmap_sigval(out_Triple, "./triple_sum_sigval.pdf")
heatmap_sigval(out_DCIS,"./dcis_sum_sigval.pdf")
heatmap_sigval(out_Invasive, "./invasive_sum_sigval.pdf")
heatmap_sigval(out_Lobular, "./lobular_sum_sigval.pdf")
heatmap_sigval(out_Papillary, "./papillary_sum_sigval.pdf")
heatmap_sigval(out_Mixed, "./mixed_sum_sigval.pdf")
heatmap_sigval(out_Grade0, "./grade0_sum_sigval.pdf")
heatmap_sigval(out_Grade1, "./grade1_sum_sigval.pdf")
heatmap_sigval(out_Grade2, "./grade2_sum_sigval.pdf")
heatmap_sigval(out_Grade3, "./grade3_sum_sigval.pdf")


# Get the list of unique classes
unique_classes <- unique(cells$class)

# Check the unique classes
print(unique_classes)

# Filter cells for "Normal" type and group by imageNb
grouped_data <- cells %>%
  filter(type == "Tumour") %>%
  group_by(imageNb) %>%
  summarise(classes_present = list(unique(class)))

# Print grouped data to check if groups are correct
print(grouped_data)

# Correctly apply filter to handle checking each row
imageNb_with_all_classes <- grouped_data %>%
  filter(sapply(classes_present, function(x) all(unique_classes %in% x))) %>%
  pull(imageNb)

# Output the imageNb that contain every class
print(imageNb_with_all_classes)



plotSpatial(spe_Normal[,spe_Normal$imageNb == "Mx_BEGIN TMA NORMAL_G13.ome_4032_4928_4544_5440.tiff"],
            node_color_by = "class", 
            draw_edges = TRUE,
            colPairName = "expansion_interaction_graph", 
            img_id = "imageNb", 
            node_size_fix = 10)+
  scale_color_manual(values = class_colors)+
  theme_minimal()

plotSpatial(spe_Normal[,spe_Normal$imageNb == "Mx_BEGIN TMA NORMAL_G13.ome_4032_4928_4544_5440.tiff"],
            node_color_by = "cn_celltypes", 
            draw_edges = TRUE,
            colPairName = "expansion_interaction_graph", 
            img_id = "imageNb", 
            node_size_fix = 10)

plotSpatial(spe_Edge[,spe_Edge$imageNb == "Mx_BEGIN TMA EDGE_I05.ome_3136_2240_3648_2752.tiff"],
            node_color_by = "class", 
            draw_edges = TRUE,
            colPairName = "expansion_interaction_graph", 
            img_id = "imageNb", 
            node_size_fix = 10)+
  scale_color_manual(values = class_colors)+
  theme_minimal()

plotSpatial(spe_Edge[,spe_Edge$imageNb == "Mx_BEGIN TMA EDGE_I05.ome_3136_2240_3648_2752.tiff"],
            node_color_by = "cn_celltypes", 
            draw_edges = TRUE,
            colPairName = "expansion_interaction_graph", 
            img_id = "imageNb", 
            node_size_fix = 10)+
  theme_minimal()



plotSpatial(spe_Tumour[,spe_Tumour$imageNb == "Mx_BEGIN TMA TUMOUR_B11.ome_4480_896_4992_1408.tiff"],
            node_color_by = "class", 
            draw_edges = TRUE,
            colPairName = "expansion_interaction_graph", 
            img_id = "imageNb", 
            node_size_fix = 10) +
  scale_color_manual(values = class_colors)+
  theme_minimal()

plotSpatial(spe_Tumour[,spe_Tumour$imageNb == "Mx_BEGIN TMA TUMOUR_B11.ome_4480_896_4992_1408.tiff"],
            node_color_by = "cn_celltypes", 
            draw_edges = TRUE,
            colPairName = "expansion_interaction_graph", 
            img_id = "imageNb", 
            node_size_fix = 10)+
  theme_minimal()

plotSpatial(spe_LuminalA[,spe_LuminalA$imageNb == "Mx_SA TMA 5_J21.ome_448_2240_960_2752.tiff"],
            node_color_by = "class", 
            draw_edges = TRUE,
            colPairName = "expansion_interaction_graph", 
            img_id = "imageNb", 
            node_size_fix = 10)

plotSpatial(spe_LuminalA[,spe_LuminalA$imageNb == "Mx_SA TMA 5_J21.ome_448_2240_960_2752.tiff"],
            node_color_by = "cn_celltypes", 
            draw_edges = TRUE,
            colPairName = "expansion_interaction_graph", 
            img_id = "imageNb", 
            node_size_fix = 10)

plotSpatial(spe_LuminalB[,spe_LuminalB$imageNb == "Mx_Southampton TMA Research Block 1_A03.ome_2688_896_3200_1408.tiff"],
            node_color_by = "class", 
            draw_edges = TRUE,
            colPairName = "expansion_interaction_graph", 
            img_id = "imageNb", 
            node_size_fix = 10)

plotSpatial(spe_LuminalB[,spe_LuminalB$imageNb == "Mx_Southampton TMA Research Block 1_A03.ome_2688_896_3200_1408.tiff"],
            node_color_by = "cn_celltypes", 
            draw_edges = TRUE,
            colPairName = "expansion_interaction_graph", 
            img_id = "imageNb", 
            node_size_fix = 10)

#Distance to proliferating cells!

spe <- minDistToCells(spe,
                      x_cells = spe$proliferating == "Positive",
                      coords = c("Pos_X","Pos_Y"),
                      img_id = 'imageNb')



print(plotSpatial(spe[,spe$imageNb == "Mx_BEGIN TMA EDGE_B11.ome"],
                  img_id = "imageNb",
                  node_color_by = "distToCells",
                  scales = "free") +
        scale_color_viridis())

counts <- t(assay(sce, "exprs"))
metadata <- as.data.frame(colData(spe))
plotting <- merge(cells,metadata,by = 'imageNb')

ggplot(plotting[plotting$imageNb== "Mx_BEGIN TMA NORMAL_B09.ome.tif",],aes(x = Pos_X, y = Pos_Y, color = class, shape = proliferating,size = area))+
  geom_point()+
  scale_color_manual(values = class_colors)+
  theme_pubr()

ggplot(plotting[plotting$imageNb== "Mx_BEGIN TMA NORMAL_B09.ome.tif",],aes(x = Pos_X, y = Pos_Y, color = cn_celltypes, shape = proliferating,size = area))+
  geom_point()+
  theme_pubr()

ggplot(plotting[plotting$imageNb== 'Mx_BEGIN TMA TUMOUR_B18.ome.tif',],aes(x = Pos_Y, y = Pos_X, color = class, shape = proliferating,size = area))+
  geom_point()+
  scale_color_manual(values = class_colors)+
  theme_pubr()


ggplot(plotting[plotting$imageNb== "Mx_BEGIN TMA NORMAL_B09.ome.tif",],aes(x = Pos_X, y = Pos_Y, color = distToCells, shape = proliferating,size = area))+
  geom_point()+
  scale_color_viridis()+
  theme_pubr()

ggplot(plotting[plotting$imageNb== "Mx_BEGIN TMA TUMOUR_B18.ome.tif",],aes(x = Pos_Y, y = Pos_X, color = distToCells, shape = proliferating,size = area))+
  geom_point()+
  scale_color_viridis()+
  theme_pubr()

grouped_plotting <- plotting %>% 
  group_by(class,cn_celltypes) %>%
  summarize('mean dist' = mean(distToCells, na.rm = TRUE))

ggplot(grouped_plotting, aes(x = Histology_description, y = `mean dist`, fill = class))+
  geom_bar(stat = 'identity',position = 'dodge')+
  scale_fill_manual(values = class_colors)

ggplot(grouped_plotting, aes(x = class, y = cn_celltypes,fill = `mean dist`))+
  geom_tile()+
  scale_fill_viridis_c()


ggplot(plotting, aes(x = cn_celltypes, y = distToCells, fill = class))+
  geom_violin()

ggplot(cells[cells$imageNb== "Mx_BEGIN TMA NORMAL_B09.ome.tif",], aes(x = Pos_X, y = Pos_Y,color = ))+
         geom_density2d_filled()




ggplot(cells[cells$imageNb == "Mx_BEGIN TMA NORMAL_B09.ome.tif",], aes(x = Pos_X, y = Pos_Y)) +
  # Density plot for positive proliferating cells
  geom_density2d_filled(data = cells[cells$imageNb == "Mx_BEGIN TMA NORMAL_B09.ome.tif" & cells$proliferating == "Positive",]) +
  # Points for all cells, colored by distToCells and sized by area
  geom_point(data = cells[cells$imageNb == "Mx_BEGIN TMA NORMAL_B09.ome.tif" & cells$proliferating != "Positive",],
             aes(color = class, size = area)) +
  # Points for proliferating cells, colored red
  geom_point(data = cells[cells$imageNb == "Mx_BEGIN TMA NORMAL_B09.ome.tif" & cells$proliferating == "Positive",],
             aes(size = area), color = "green",shape = 18) +
  # Density lines for positive proliferating cells
  geom_density_2d(data = cells[cells$imageNb == "Mx_BEGIN TMA NORMAL_B09.ome.tif" & cells$proliferating == "Positive",]) +
  scale_fill_viridis(option = "B",discrete = TRUE) +
  scale_color_manual(values = class_colors)+
  theme_pubr() +
  labs(x = "X Position",y = "Y Position",color = "Cell Type",size = "Area",fill = "Density")


ggplot(cells[cells$imageNb == "Mx_BEGIN TMA NORMAL_B09.ome.tif",], aes(x = Pos_X, y = Pos_Y)) +
  # This would create filled polygons that we use for create our polygon
  geom_density2d_filled() +
  geom_point(data = cells[cells$imageNb == "Mx_BEGIN TMA NORMAL_B09.ome.tif" & cells$proliferating != "Positive",],
             aes(color = class, size = area)) +
  # Points for proliferating cells, colored red
  geom_point(data = cells[cells$imageNb == "Mx_BEGIN TMA NORMAL_B09.ome.tif" & cells$proliferating == "Positive",],
             aes(size = area), color = "green",shape = 18) +
  geom_density_2d()+
  scale_fill_viridis(option = "H",discrete = TRUE)+
  scale_color_manual(values = class_colors)+
  theme_pubr() +
  labs(x = "X Position",y = "Y Position",color = "Cell Type",size = "Area",fill = "Density")



#-----------------------------------------------------------#
ggplot(cells[cells$imageNb == "Mx_BEGIN TMA EDGE_F09.ome.tif",], aes(x = Pos_X, y = Pos_Y)) +
  # Density plot for positive proliferating cells
  geom_density2d_filled(data = cells[cells$imageNb == "Mx_BEGIN TMA EDGE_F09.ome.tif" & cells$proliferating == "Positive",]) +
  # Points for all cells, colored by distToCells and sized by area
  geom_point(data = cells[cells$imageNb == "Mx_BEGIN TMA EDGE_F09.ome.tif" & cells$proliferating != "Positive",],
             aes(color = ANAX1, size = area)) +
  # Points for proliferating cells, colored red
  geom_point(data = cells[cells$imageNb == "Mx_BEGIN TMA EDGE_F09.ome.tif" & cells$proliferating == "Positive",],
             aes(size = area,color = ANAX1),shape = 18) +
  # Density lines for positive proliferating cells
  geom_density_2d(data = cells[cells$imageNb == "Mx_BEGIN TMA EDGE_F09.ome.tif" & cells$proliferating == "Positive",]) +
  scale_fill_viridis(option = "B",discrete = TRUE) +
  
  theme_pubr() +
  labs(x = "X Position",y = "Y Position",color = "Cell Type",size = "Area",fill = "Density")


merged <- cbind(counts, plotting)

ggplot(merged,aes(x = distToCells, y = ANAX1))+
  geom_point()+
  theme_minimal() +
  coord_cartesian() +
  labs(x = "ANAX1", y = "Distance to Cells")


ggplot(merged,aes(x = ANAX1))+
  geom_histogram()



normal <- merged[merged$Type == "Edge",]

plot(normal$distToCells,normal$ANAX1)

# Perform kernel density estimation
kde <- kde2d(merged$Pos_X, merged$Pos_Y, n = 100)
density_df <- data.frame(expand.grid(x = kde$x, y = kde$y))
density_df$z <- c(kde$z)


# Merge density values with original data
plotting_with_density <- merge(plotting, density_df, by.x = c("Pos_X", "Pos_Y"), by.y = c("x", "y"), all.x = TRUE)

# Create density bands
plotting_with_density <- plotting_with_density %>%
  mutate(density_band = cut(z, breaks = quantile(z, probs = 0:5 / 5), include.lowest = TRUE))

# Calculate the mean of 'area' for each density band
mean_by_band <- plotting_with_density %>%
  group_by(density_band) %>%
  summarise(mean_area = mean(area, na.rm = TRUE), .groups = 'drop')

print(mean_by_band)

# Plot to visualize (optional)
ggplot(plotting_with_density, aes(x = Pos_X, y = Pos_Y, fill = density_band)) +
  geom_tile() +
  geom_point(aes(size = area, color = z)) +
  scale_fill_viridis_c() +
  scale_color_viridis_c() +
  theme_minimal()