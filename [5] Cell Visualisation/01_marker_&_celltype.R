library(ggplot2)
library(ggraph)
library(viridis)
library(tidyverse)
library(pheatmap)
library(scales)
library(BiocParallel)
library(data.table)
library(RColorBrewer)
library(dplyr)
library(pals)
library(janitor)
library(ggpubr)
library(rstatix)
library(ggbeeswarm)
library(glue)
library(dittoSeq)
library(scater)
library(patchwork)
library(cowplot)
library(gplots)

set.seed(123)
class_colors <- c(Luminal = "#51abcb", BnL = "#f8766d",Basal = "red", `Blood Vessel` = "#a3a500", LnB = "#00b0f6",Undefined = "grey")
Type <- c("#3274A1","#FF812C","#3A923A")
Diag <- c("#A9DE00","#60C943","#3EA074","#006769")
Hist <- c("#FFBEFF", "#FF70FF", "#CB64B8", "#C964D1", "#A05EA3")
Grade <- c("#1E90FF", "#4682B4", "#5F9EA0", "#B0C4DE")
cells <- read.csv("E:/FINAL_SLICE_IMG_MERGED_PATIENTDATA.csv")

cells$Type <- factor(cells$Type, levels = c("Normal", "Edge", "Tumour"))
cells$diagnosis <- factor(cells$diagnosis, levels = c("Luminal A","Luminal B","HER2 enriched","Triple negative"))
cells$Histology_description <- factor(cells$Histology_description, levels = c("DCIS","Invasive ductal","Lobular","Papillary","Mixed"))
cells$GRADE <- factor(cells$GRADE)
cells$class <- factor(cells$class, levels = c("Basal","BnL","Luminal","LnB","Blood Vessel","Undefined"))

#----------------------------------------------------------------------------------------------------------#
#                                        Marker expression Plots                                           #
#----------------------------------------------------------------------------------------------------------#

# Ensure the necessary directories exist or create them
create_dir_if_not_exists <- function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}


create_violin_plot <- function(marker, cell_type) {
  # setwd("C:/Users/miles/GitHub/dissertation//Figure 4(Markers)")
  dir <- file.path(cell_type, marker)
  create_dir_if_not_exists(dir)
  
  
  plt1 <- ggplot(cells[(cells$class == cell_type & !(is.na(cells$Type))),],aes(x = Type, y = .data[[marker]],fill = Type))+
    geom_violin()+
    geom_boxplot(width=.1, outlier.shape=NA,alpha = .2) +
    scale_fill_manual(values = Type)+
    labs(x = "Classification", y = paste(marker," marker expression(a.u.)")) +
    scale_y_continuous(breaks = seq(0, max(cells[,marker]), by = (max(cells[,marker])) / 10)) +
    theme_pubr()+
    stat_compare_means(method = "wilcox",label = "p.signif", comparisons = list(c("Normal","Edge"),c("Normal","Tumour"),c("Edge","Tumour")))+
    theme(axis.text = element_text(face="bold", size = 15),axis.title.y = element_text(face="bold", size = 15),axis.title.x = element_text(face="bold", size = 15),legend.title = element_text(face = "bold"))
  
  
  ggsave(filename = file.path(dir, paste(marker, "Type_Classification.pdf", sep = "_")), plot = plt1, width = 10, height = 7)
  
  plt2 <-ggplot(cells[(cells$class == cell_type),],aes(x = diagnosis, y = .data[[marker]],fill = diagnosis))+
    geom_violin()+
    geom_boxplot(width=.1, outlier.shape=NA,alpha = .2) +
    labs(x = "Patient Diagnosis", y = paste(marker," marker expression(a.u.)")) +
    scale_fill_manual(name = "Diagnosis",values = Diag)+
    scale_y_continuous(breaks = seq(0, max(cells[,marker]), by = (max(cells[,marker]) / 10))) +
    theme_pubr()+
    stat_compare_means(method = "kruskal.test",aes(label = paste0("p",scales::label_pvalue()(..p..),"****")),label.x = 0.7, label.y = 11)+
    theme(axis.text = element_text(face="bold", size = 15),axis.title.y = element_text(face="bold", size = 15),axis.title.x = element_text(face="bold", size = 15),legend.title = element_text(face = "bold"))
  
  ggsave(filename = file.path(dir, paste(marker, "Patient_Diagnosis.pdf", sep = "_")), plot = plt2, width = 10, height = 7)
  
  plt3 <-ggplot(cells[(cells$class == cell_type),],aes(x = Histology_description, y = .data[[marker]],fill = Histology_description))+
    geom_violin()+
    geom_boxplot(width=.1, outlier.shape=NA,alpha = .2) +
    labs(x = "Tumour Histology", y = paste(marker," marker expression(a.u.)")) +
    scale_fill_manual(name = "Histology",values = Hist)+
    scale_y_continuous(breaks = seq(0, max(cells[,marker]), by = (max(cells[,marker])) / 10)) +
    theme_pubr()+
    stat_compare_means(method = "kruskal.test",aes(label = paste0("p",scales::label_pvalue()(..p..),"****")),label.x = 0.7, label.y = 11)+
    theme(axis.text = element_text(face="bold", size = 15),axis.title.y = element_text(face="bold", size = 15),axis.title.x = element_text(face="bold", size = 15),legend.title = element_text(face = "bold"))
  
  ggsave(filename = file.path(dir, paste(marker, "Tumour_Histology.pdf", sep = "_")), plot = plt3, width = 10, height = 7)
  
  plt4 <-ggplot(cells[(cells$class == cell_type),],aes(x = GRADE, y = .data[[marker]],fill = GRADE))+
    geom_violin()+
    geom_boxplot(width=.1, outlier.shape=NA,alpha = .2) +
    labs(x = "Tumour Grade", y = paste(marker," marker expression(a.u.)")) +
    scale_fill_manual(name = "Grade",values = Grade)+
    scale_y_continuous(breaks = seq(0, max(cells[,marker]), by = (max(cells[,marker])) / 10)) +
    theme_pubr()+
    stat_compare_means(method = "kruskal.test",aes(label = paste0("p",scales::label_pvalue()(..p..),"***")),label.x = 0.7, label.y = 11)+
    theme(axis.text = element_text(face="bold", size = 15),axis.title.y = element_text(face="bold", size = 15),axis.title.x = element_text(face="bold", size = 15),legend.title = element_text(face = "bold"))
  
  ggsave(filename = file.path(dir, paste(marker, "Tumour_Grade.pdf", sep = "_")), plot = plt4, width = 10, height = 7)
  
  print(plt1)
  print(plt2)
  print(plt3)
  print(plt4)
}

cells[cells$Full_Filename == "Mx_BEGIN TMA TUMOUR_J17.ome.tif",]

cells[cells$Full_Filename == "Mx_BEGIN TMA TUMOUR_J17.ome.tif",]


create_violin_plot("area","Luminal")
create_violin_plot("orientation","Luminal")

create_violin_plot("ANAX1","LnB")
create_violin_plot("ANAX1","BnL")
create_violin_plot("ANAX1","Basal")

for (marker in markers_columns) {
  create_violin_plot(marker,"Luminal")
}


for (marker in markers_columns) {
  create_violin_plot(marker,"BnL")
}


for (marker in markers_columns) {
  create_violin_plot(marker,"LnB")
}


for (marker in markers_columns) {
  create_violin_plot(marker,"Basal")
}

for (marker in markers_columns) {
  create_violin_plot(marker,"Undefined")
}

for (marker in markers_columns) {
  create_violin_plot(marker,"Blood Vessel")
}

#----------------------------------------------------------------------------------------------------------#
#                                        Cell Count Plots                                                  #
#----------------------------------------------------------------------------------------------------------#


cell_counts <- table(cells$class)#
# Generate the table
cell_counts <- as.data.frame(table(cells$class))

# Rename columns for clarity
colnames(cell_counts) <- c("class", "count")

ggplot(cell_counts, aes(x = class,  y = count,fill = class)) +
  geom_bar(stat = "identity",position = "dodge",color = "black") +
  geom_text(aes(label = scales::comma(count)), vjust = -0.5, size = 3.5) + 
  labs(x = "Cell Type", y = "Cell Count", fill = "Cell type") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = scales::comma,name="Cell Counts",limits = c(0, max(cell_counts$count) * 1.1))+
  scale_fill_manual(values = class_colors)+
  ggtitle("Proportion of Cells in Each Image")+
  theme_pubr()

# 
# 
# ggplot(cells, aes(x = Type, fill = class))+
#   geom_bar(stat = "count",position = "dodge")
# 
# 
# cell_counts <- cell_counts %>%
#   group_by(diagnosis) %>%
#   mutate(Percentage = Count / sum(Count))
# 
# cell_counts$Percentage <- as.numeric(cell_counts$Percentage)
# 
# 
# ggplot(cell_counts, aes(x = diagnosis, y = Percentage, fill = class)) +
#   geom_bar(stat = "identity") +
#   scale_y_continuous(labels = scales::percent_format(scale = 100),breaks = seq(0, 1, by = 0.1)) +
#   labs(x = "Image", y = "Percentage of Cell Types", fill = "Cell Class") +
#   scale_fill_manual(fill())
#   theme_pubr() +
#   ggtitle("Percentage proportion of cell types")
# 
# 
# ggplot(cell_counts, aes(x = Type, y = Percentage, fill = class)) +
#   geom_bar(stat = "identity") +
#   scale_y_continuous(labels = scales::percent_format(scale = 1)) + # Correct percentage scaling
#   labs(x = "Diagnosis", y = "Percentage", fill = "Cell Type") +
#   theme_pubr() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   ggtitle("Proportion of Cells in Each Diagnosis") +
#   facet_wrap(~ class, scales = "free_y")


diagnosis_summary <- cells %>%
  group_by(diagnosis, class) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(diagnosis) %>%
  mutate(percentage = (count / sum(count)) * 100) %>%
  ungroup() %>%
  select(diagnosis, class, percentage)

# Plot for Diagnosis
ggplot(diagnosis_summary, aes(x = diagnosis, y = percentage, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)),
            position = position_dodge(width = 0.9),
            vjust = -0.25, size = 2) +
  coord_cartesian(ylim = c(0, 100))+
  labs(title = "Percentage of Cell Types by Diagnosis",
       x = "Diagnosis",
       y = "Percentage",
       fill = "Class") +
  scale_y_continuous(name="Percentage of cells",breaks = seq(0, 100, by = 10) ,labels = function(x) paste0(x,"%"))+
  scale_fill_manual(values = class_colors)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_pubr()

histology_summary <- cells %>%
  group_by(Histology_description, class) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(Histology_description) %>%
  mutate(percentage = (count / sum(count)) * 100) %>%
  ungroup() %>%
  select(Histology_description, class, percentage)

# Plot for Histology
ggplot(histology_summary, aes(x = Histology_description, y = percentage, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)),
            position = position_dodge(width = 0.9),
            vjust = -0.25, size = 2) +
  coord_cartesian(ylim = c(0, 100))+
  scale_fill_manual(values = class_colors)+
  scale_y_continuous(name="Percentage of cells",breaks = seq(0, 100, by = 10) ,labels = function(x) paste0(x,"%"))+
  labs(title = "Percentage of Cell Types by Histology",
       x = "Histology",
       y = "Percentage",
       fill = "Class") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_pubr()

# Generate the summary table for type
type_summary <- cells[!(is.na(cells$Type)),] %>%
  group_by(Type, class) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(Type) %>%
  mutate(percentage = (count / sum(count)) * 100) %>%
  ungroup() %>%
  select(Type, class, percentage)

# Plot for Type
ggplot(type_summary, aes(x = Type, y = percentage, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)),
            position = position_dodge(width = 0.9),
            vjust = -0.25, size = 3) +
  coord_cartesian(ylim = c(0, 100))+
  scale_fill_manual(values = class_colors)+
  scale_y_continuous(name="Percentage of cells",breaks = seq(0, 100, by = 10) ,labels = function(x) paste0(x,"%"))+
  labs(title = "Percentage of Cell Types by Type",
       x = "Type",
       y = "Percentage",
       fill = "Class") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_pubr()




# Function to summarize data with raw counts
summarize_counts <- function(df, group_col, class_col) {
  df %>%
    group_by(across(all_of(group_col)), class) %>%
    summarise(count = n(), .groups = 'drop') %>%
    ungroup()
}

# Function to plot data with raw counts and save the plot
plot_data_counts <- function(summary_df, x_col, title, p_value, file_name) {
  plot <- summary_df %>%
    group_by(across(all_of(x_col))) %>%
    mutate(percentage = (count / sum(count)) * 100) %>%
    ungroup() %>%
    ggplot(aes_string(x = x_col, y = "percentage", fill = "class")) +
    geom_bar(stat = "identity", position = "dodge",color = "black") +
    geom_text(aes(label = sprintf("%.1f%%", percentage)), 
              position = position_dodge(width = 0.9), 
              vjust = -0.25, size = 2, fontface = "bold") +
    coord_cartesian(ylim = c(0, 100)) +
    labs(title = title,
         x = x_col,
         y = "Percentage",
         fill = "Class") +
    scale_y_continuous(name = "Percentage of cells", breaks = seq(0, 100, by = 10), labels = function(x) paste0(x, "%")) +
    scale_fill_manual(values = class_colors) +
    annotate("text", x = 1, y = 90, size = 3, label = insight::format_p(p_value, stars = TRUE)) +
    theme_pubr() +
    theme(axis.text = element_text(face="bold", size = 15),
          axis.title.y = element_text(face="bold", size = 15),
          axis.title.x = element_text(face="bold", size = 15),
          legend.title = element_text(face = "bold"))
  
  ggsave(filename = file_name, plot = plot, width = 10, height = 7)
}

# Summarize and plot data for Diagnosis
diagnosis_counts <- summarize_counts(cells, "diagnosis", "class")
diagnosis_counts_wide <- diagnosis_counts %>% pivot_wider(names_from = class, values_from = count)
chi_square_result_diagnosis <- chisq.test(diagnosis_counts_wide[, -1])
plot_data_counts(diagnosis_counts, "diagnosis", "Percentage of Cell Types by Diagnosis", chi_square_result_diagnosis$p.value, "diagnosis_plot.pdf")

# Summarize and plot data for Histology
histology_counts <- summarize_counts(cells, "Histology_description", "class")
histology_counts_wide <- histology_counts %>% pivot_wider(names_from = class, values_from = count)
chi_square_result_histology <- chisq.test(histology_counts_wide[, -1])
plot_data_counts(histology_counts, "Histology_description", "Percentage of Cell Types by Histology", chi_square_result_histology$p.value, "histology_plot.pdf")

# Summarize and plot data for Type
type_counts <- cells %>% filter(!is.na(Type)) %>% summarize_counts("Type", "class")
type_counts_wide <- type_counts %>% pivot_wider(names_from = class, values_from = count)
chi_square_result_type <- chisq.test(type_counts_wide[, -1])
plot_data_counts(type_counts, "Type", "Percentage of Cell Types by Type", chi_square_result_type$p.value, "type_plot.pdf")

# Summarize and plot data for Grade
grade_counts <- summarize_counts(cells, "GRADE", "class")
grade_counts_wide <- grade_counts %>% pivot_wider(names_from = class, values_from = count)
chi_square_result_grade <- chisq.test(grade_counts_wide[, -1])
plot_data_counts(grade_counts, "GRADE", "Percentage of Cell Types by Grade", chi_square_result_grade$p.value, "grade_plot.pdf")

#----------------------------------------------------------------------------------------------------------#
#                                             Proliferating                                               #
#----------------------------------------------------------------------------------------------------------#

cell_counts <- table(cells$class[cells$proliferating == "Positive"])
# Generate the table
cell_counts <- as.data.frame((cell_counts))

# Rename columns for clarity
colnames(cell_counts) <- c("class", "count")

ggplot(cell_counts, aes(x = class,  y = count,fill = class)) +
  geom_bar(stat = "identity",position = "dodge",color = "black") +
  geom_text(aes(label = scales::comma(count)), vjust = -0.5, size = 3.5) + 
  labs(x = "Cell Type", y = "Cell Count", fill = "Cell type") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = scales::comma,name="Cell Counts",limits = c(0, max(cell_counts$count) * 1.1))+
  scale_fill_manual(values = class_colors)+
  ggtitle("Proportion of proliferating cells")+
  theme_pubr()


# Proliferating Cells by Type
proliferating_type_summary <- cells[!(is.na(cells$type)),] %>%
  group_by(Type, proliferating, class) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(Type) %>%
  mutate(percentage = (count / sum(count)) * 100) %>%
  ungroup()

type_table <- table(cells$type, cells$proliferating)
type_chisq <- chisq.test(type_table)
type_p_value <- type_chisq$p.value

ggplot(proliferating_type_summary[proliferating_type_summary$proliferating == "Positive",], 
       aes(x = Type, y = percentage, fill = class)) +
  geom_bar(stat = "identity", position = "dodge",color = "black") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.25, size = 3) +
  coord_cartesian(ylim = c(0, 10)) +
  labs(title = "Percentage of Proliferating Cells by Type",
       x = "Type",
       y = "Percentage",
       fill = "Class") +
  scale_fill_manual(values = class_colors) +
  scale_y_continuous(name = "Percentage of Proliferating Cells",
                     breaks = seq(0, 10, by = 1), 
                     labels = function(x) paste0(x, "%")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  annotate("text", x = 1, y = 9, size = 3, 
           label = insight::format_p(type_p_value, stars = TRUE))+
  theme_pubr() +
  theme(axis.text = element_text(face="bold", size = 15),
        axis.title.y = element_text(face="bold", size = 15),
        axis.title.x = element_text(face="bold", size = 15),
        legend.title = element_text(face = "bold"))


# Proliferating Cells by Diagnosis
proliferating_summary <- cells[!(is.na(cells$diagnosis)),] %>%
  group_by(diagnosis, proliferating, class) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(diagnosis) %>%
  mutate(percentage = (count / sum(count)) * 100) %>%
  ungroup()

diagnosis_table <- table(cells$diagnosis, cells$proliferating)
diagnosis_chisq <- chisq.test(diagnosis_table)
diagnosis_p_value <- diagnosis_chisq$p.value

ggplot(proliferating_summary[proliferating_summary$proliferating == "Positive",], 
       aes(x = diagnosis, y = percentage, fill = class)) +
  geom_bar(stat = "identity", position = "dodge",color = "black") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.25, size = 3) +
  coord_cartesian(ylim = c(0, 10)) +
  labs(title = "Percentage of Proliferating Cells by Diagnosis",
       x = "Diagnosis",
       y = "Percentage",
       fill = "Class") +
  scale_fill_manual(values = class_colors) +
  scale_y_continuous(name = "Percentage of Proliferating Cells",
                     breaks = seq(0, 10, by = 1), 
                     labels = function(x) paste0(x, "%")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  annotate("text", x = 1, y = 9, size = 3, 
           label = insight::format_p(diagnosis_p_value, stars = TRUE))+
  theme_pubr() +
  theme(axis.text = element_text(face="bold", size = 15),
        axis.title.y = element_text(face="bold", size = 15),
        axis.title.x = element_text(face="bold", size = 15),
        legend.title = element_text(face = "bold"))

# Proliferating Cells by Histology
proliferating_histology_summary <- cells[!(is.na(cells$Histology_description)),] %>%
  group_by(Histology_description, proliferating, class) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(Histology_description) %>%
  mutate(percentage = (count / sum(count)) * 100) %>%
  ungroup()

histology_table <- table(cells$Histology_description, cells$proliferating)
histology_chisq <- chisq.test(histology_table)
histology_p_value <- histology_chisq$p.value

ggplot(proliferating_histology_summary[proliferating_histology_summary$proliferating == "Positive",], 
       aes(x = Histology_description, y = percentage, fill = class)) +
  geom_bar(stat = "identity", position = "dodge",color = "black") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.25, size = 3) +
  coord_cartesian(ylim = c(0, 10)) +
  labs(title = "Percentage of Proliferating Cells by Histology",
       x = "Histology",
       y = "Percentage",
       fill = "Class") +
  scale_fill_manual(values = class_colors) +
  scale_y_continuous(name = "Percentage of Proliferating Cells",
                     breaks = seq(0, 10, by = 1), 
                     labels = function(x) paste0(x, "%")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  annotate("text", x = 1, y = 9, size = 3, 
           label = insight::format_p(histology_p_value, stars = TRUE))+
  theme_pubr() +
  theme(axis.text = element_text(face="bold", size = 15),
        axis.title.y = element_text(face="bold", size = 15),
        axis.title.x = element_text(face="bold", size = 15),
        legend.title = element_text(face = "bold"))

# Proliferating Cells by Grade
proliferating_grade_summary <- cells[!(is.na(cells$GRADE)),] %>%
  group_by(GRADE, proliferating, class) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(GRADE) %>%
  mutate(percentage = (count / sum(count)) * 100) %>%
  ungroup()

grade_table <- table(cells$GRADE, cells$proliferating)
grade_chisq <- chisq.test(grade_table)
grade_p_value <- grade_chisq$p.value

ggplot(proliferating_grade_summary[proliferating_grade_summary$proliferating == "Positive",], 
       aes(x = GRADE, y = percentage, fill = class)) +
  geom_bar(stat = "identity", position = "dodge",color = "black") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.25, size = 3) +
  coord_cartesian(ylim = c(0, 10)) +
  labs(title = "Percentage of Proliferating Cells by Grade",
       x = "Grade",
       y = "Percentage",
       fill = "Class") +
  scale_fill_manual(values = class_colors) +
  scale_y_continuous(name = "Percentage of Proliferating Cells",
                     breaks = seq(0, 10, by = 1), 
                     labels = function(x) paste0(x, "%")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  annotate("text", x = 1, y = 9, size = 3, 
           label = insight::format_p(grade_p_value, stars = TRUE))+
  theme_pubr() +
  theme(axis.text = element_text(face="bold", size = 15),
        axis.title.y = element_text(face="bold", size = 15),
        axis.title.x = element_text(face="bold", size = 15),
        legend.title = element_text(face = "bold"))