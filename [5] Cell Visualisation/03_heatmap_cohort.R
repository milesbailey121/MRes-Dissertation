library(tidyverse)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(pheatmap)
library(viridis)
library(dplyr)

set.seed(123)
data <- read_csv("E:/Final Data/quantification/mean_core_intensity.csv")
markers_columns <- c("CD31","SMA","Ki67","CK8","ANAX1","ECAD")
class_colors <- c(Luminal = "#51abcb", BnL = "#f8766d",Basal = "red", `Blood Vessel` = "#a3a500", LnB = "#00b0f6",Undefined = "grey")
Type <- c("#3274A1","#FF812C","#3A923A","lightgrey")
Diag <- c("#A9DE00","#60C943","#3EA074","#006769")
Hist <- c("#FFBEFF", "#FF70FF", "#CB64B8", "#C964D1", "#A05EA3")
Grade <- c("#1E90FF", "#4682B4", "#5F9EA0", "#B0C4DE")

data$type <- factor(data$type, levels = c("Normal", "Edge", "Tumour"))
data$diagnosis <- factor(data$diagnosis, levels = c("Luminal A","Luminal B","HER2 enriched","Triple negative"))
data$Histology_description <- factor(data$Histology_description, levels = c("DCIS","Invasive ductal","Lobular","Papillary","Mixed"))
data$GRADE <- factor(data$GRADE)



heatmap_data <- data[, markers_columns]
rownames(heatmap_data) <- data$imageNb

mat_col <- data.frame(group = data$type)

mat_colors <- list(group = brewer.pal(4, "Set1"))
names(mat_colors$group) <- unique(data$type)


pheatmap(
  mat               = heatmap_data,
  color             = inferno(10),
  border_color      = NA,
  show_colnames     = TRUE,
  show_rownames     = FALSE,
  cluster_rows = FALSE,
  annotation_row    = mat_col,
  annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 14,
  main              = "Default Heatmap"
)





# Extract the relevant data
heatmap_data <- data[, markers_columns]
rownames(heatmap_data) <- data$imageNb

# Create the annotation DataFrame
mat_col <- data.frame(
  Type = data$type,
  Diag = data$diagnosis,
  Hist = data$Histology_description,
  Grade = data$GRADE
)

mat_col$Type <- factor(mat_col$Type, levels = c("Normal", "Edge", "Tumour"))
mat_col$Diag <- factor(mat_col$Diag , levels = c("Luminal A","Luminal B","HER2 enriched","Triple negative"))
mat_col$Hist <- factor(mat_col$Hist, levels = c("DCIS","Invasive ductal","Lobular","Papillary","Mixed"))
mat_col$Grade <- factor(mat_col$Grade)


rownames(mat_col) <- data$imageNb

# Ensure that NA values are represented
mat_col$Type[is.na(mat_col$Type)] <- "Unclassified"  # Assign a placeholder for NA values
mat_col$Type <- factor(mat_col$Type, levels = c("Normal", "Edge", "Tumour", "Unclassified"))

# Define colors for each group, including NA
mat_colors <- list(
  Type = c("Normal" = "#3274A1", "Edge" = "#FF812C", "Tumour" = "#3A923A", "Unclassified" = "#D3D3D3"),
  Diag = c("Luminal A" = "#A9DE00", "Luminal B" = "#60C943", "HER2 enriched" = "#3EA074", "Triple negative" = "#006769"),
  Hist = c("DCIS" = "#FFBEFF", "Invasive ductal" = "#FF70FF", "Lobular" = "#CB64B8", "Papillary" = "#C964D1", "Mixed" = "#A05EA3"),
  Grade = c("0" = "#1E90FF", "1" = "#4682B4", "2" = "#5F9EA0", "3" = "#B0C4DE")
)

# Map colors to the actual levels in the data
names(mat_colors$Type) <- unique(data$type)
names(mat_colors$Diag) <- unique(data$diagnosis)
names(mat_colors$Hist) <- unique(data$Histology_description)
names(mat_colors$Grade) <- unique(data$GRADE)



# Convert categorical annotations to numeric values
annotation_numeric <- data.frame(
  Type = as.numeric(factor(data$type)),
  Diag = as.numeric(factor(data$diagnosis)),
  Hist = as.numeric(factor(data$Histology_description)),
  Grade = as.numeric(factor(data$GRADE))
)
rownames(annotation_numeric) <- data$imageNb

# Calculate the distance matrix based on annotations
annotation_dist <- dist(annotation_numeric)

# Perform hierarchical clustering based on the annotation distance
annotation_clustering <- hclust(annotation_dist)

# Plot the heatmap with clustering based on annotations
pheatmap(
  mat               = heatmap_data,
  color             = inferno(10),           # Color palette for heatmap
  border_color      = NA,                    # No borders around cells
  show_colnames     = TRUE,                  # Show column names (markers)
  show_rownames     = FALSE,                 # Hide row names (filenames)
  annotation_row    = mat_col,               # Add row annotations
  annotation_colors = mat_colors,            # Add colors for annotations
  drop_levels       = TRUE,                  # Drop unused levels
  fontsize          = 14,                    # Font size for labels
  cluster_rows      = annotation_clustering, # Apply custom row clustering
  cluster_cols      = TRUE,                 # Disable column clustering
  main              = "Heatmap Clustered by Annotations"  # Title of the heatmap
)

# Plot the heatmap with annotations
pheatmap(
  mat               = heatmap_data,
  color             = inferno(10),           # Color palette for heatmap
  border_color      = NA,                    # No borders around cells
  show_colnames     = TRUE,                  # Show column names (markers)
  show_rownames     = FALSE,
  # cluster_rows = FALSE,# Hide row names (filenames)
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_row    = mat_col,               # Add row annotations
  annotation_colors = mat_colors,            # Add colors for annotations
  drop_levels       = TRUE,                  # Drop unused levels
  fontsize          = 14,                    # Font size for labels
  main              = "Grouped Heatmap"      # Title of the heatmap
)


ggplot(data[!(is.na(data$type)),], aes(x = type, y = ANAX1, fill = type)) +
  geom_violin(trim = FALSE,alpha = 0.5) +
  geom_beeswarm(shape = 16, color = "black", alpha = 1) +
  geom_boxplot(width=.1, outlier.shape=NA,alpha = .2) +
  scale_y_continuous(breaks = ) +
  stat_compare_means(method = "wilcox",label = "p.signif", comparisons = list(c("Normal","Edge"),c("Normal","Tumour"),c("Edge","Tumour"))) +
  scale_fill_manual(values = Type) +
  labs(title = "Violin Plot of ANAX1 vs Type", x = "Type", y = "ANAX1 Intensity") +
  theme_pubr()

ggplot(data, aes(x = diagnosis, y = ANAX1, fill = diagnosis)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_beeswarm(shape = 16, color = "black", alpha = 1) +
  geom_boxplot(width = .1, outlier.shape = NA, alpha = .2) +
  scale_y_continuous(breaks = ) +
  stat_compare_means(method = "wilcox", label = "p.signif", comparisons = list(
    c("Luminal A", "Luminal B"),
    c("Luminal A", "HER2 enriched"),
    c("Luminal A", "Triple negative"),
    c("Luminal B", "HER2 enriched"),
    c("Luminal B", "Triple negative"),
    c("HER2 enriched", "Triple negative")
  )) +
  scale_fill_manual(values = Diag) +
  labs(title = "Violin Plot of ANAX1 vs Diagnosis", x = "Diagnosis", y = "ANAX1 Intensity") +
  theme_pubr()


ggplot(data, aes(x = Histology_description, y = ANAX1, fill = Histology_description)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_beeswarm(shape = 16, color = "black", alpha = 1) +
  geom_boxplot(width = .1, outlier.shape = NA, alpha = .2) +
  scale_y_continuous(breaks = ) +
  stat_compare_means(method = "wilcox", label = "p.signif", comparisons = list(
    c("DCIS", "Invasive ductal"),
    c("DCIS", "Lobular"),
    c("DCIS", "Papillary"),
    c("DCIS", "Mixed"),
    c("Invasive ductal", "Lobular"),
    c("Invasive ductal", "Papillary"),
    c("Invasive ductal", "Mixed"),
    c("Lobular", "Papillary"),
    c("Lobular", "Mixed"),
    c("Papillary", "Mixed")
  )) +
  scale_fill_manual(values = Hist) +
  labs(title = "Violin Plot of ANAX1 vs Histology", x = "Histology", y = "ANAX1 Intensity") +
  theme_pubr()


ggplot(data, aes(x = GRADE, y = ANAX1, fill = GRADE)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_beeswarm(shape = 16, color = "black", alpha = 1) +
  geom_boxplot(width = .1, outlier.shape = NA, alpha = .2) +
  scale_y_continuous(breaks = ) +
  stat_compare_means(method = "t.test", label = "p.signif", comparisons = list(
    c("0", "1"),
    c("0", "2"),
    c("0", "3"),
    c("1", "2"),
    c("1", "3"),
    c("2", "3")
  )) +
  scale_fill_manual(values = Grade) +
  labs(title = "Violin Plot of ANAX1 vs Grade", x = "Grade", y = "ANAX1 Intensity") +
  theme_pubr()




Type2 <- c("#3274A1","#FF812C","#3A923A","lightgrey")
Diag2 <- c("#3274A1","#A9DE00","#60C943","#3EA074","#006769")
Hist2 <- c("#3274A1","#FFBEFF", "#FF70FF", "#CB64B8", "#C964D1", "#A05EA3")
Grade2 <- c("#3274A1","#1E90FF", "#4682B4", "#5F9EA0", "#B0C4DE")



# Filter the "Normal" and "Edge" data from the original dataset
data_normal <- data %>% filter(type %in% c("Normal"))

data_normal$diagnosis <- "Normal"
data_normal$Histology_description <- "Normal" 

combined_data <- bind_rows(data_normal,data)

# Now plot with ggplot2 for each annotation, e.g., 'Diagnosis'

combined_data$diagnosis <- factor(combined_data$diagnosis, levels = c("Normal","Luminal A","Luminal B","HER2 enriched","Triple negative"))
combined_data$Histology_description <- factor(combined_data$Histology_description, levels = c("Normal","DCIS","Invasive ductal","Lobular","Papillary","Mixed"))


# Violin Plot for ANAX1 vs Diagnosis with Normal and Edge
ggplot(combined_data, aes(x = diagnosis, y = ANAX1, fill = diagnosis)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_beeswarm(shape = 16, color = "black",alpha = 1,size = 2) +
  geom_boxplot(width = .1, outlier.shape = NA, alpha = .2) +
  ylim(-2,10)+
  stat_compare_means(method = "wilcox", label = "p.signif", comparisons = list(c("Normal", "Luminal A"), c("Normal", "Luminal B"), c("Normal", "HER2 enriched"), c("Normal", "Triple negative"))) +
  scale_fill_manual(values = Diag2) +
  labs(title = "Violin Plot of ANAX1 vs Diagnosis with Normal and Edge", x = "Diagnosis", y = "ANAX1 Intensity") +
  theme_pubr()


# Plot for Histology
ggplot(combined_data, aes(x = Histology_description, y = ANAX1, fill = Histology_description)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_beeswarm(shape = 16, color = "black", alpha = 1) +
  geom_boxplot(width = .1, outlier.shape = NA, alpha = .2) +
  ylim(-2,10)+
  stat_compare_means(method = "wilcox", label = "p.signif", comparisons = list(c("Normal", "DCIS"),c("Normal", "Invasive ductal"), c("Normal", "Lobular"), c("Normal", "Papillary"), c("Normal", "Mixed"))) +
  scale_fill_manual(values = Hist2) +
  labs(title = "Violin Plot of ANAX1 vs Histology", x = "Histology Description", y = "ANAX1 Intensity") +
  theme_pubr()

# Plot for Grade
ggplot(data, aes(x = GRADE, y = ANAX1, fill = GRADE)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_beeswarm(shape = 16, color = "black", alpha = 1) +
  geom_boxplot(width = .1, outlier.shape = NA, alpha = .2) +
  ylim(-2,10)+
  scale_fill_manual(values = Grade2)+
  stat_compare_means(method = "wilcox", label = "p.signif", comparisons = list(c("0", "1"),c("0", "2"), c("0", "3"))) +
  labs(title = "Violin Plot of ANAX1 vs Grade", x = "Grade", y = "ANAX1 Intensity") +
  theme_pubr()



