
#----------------------------------------------Read Data In -----------------------------------------------#

patient_data <- read_csv("./path/to/data")

cells <- read.csv("./path/to/data")

#----------------------------------------------------------------------------------------------------------#
#                                         Patient Data                                                     #
#----------------------------------------------------------------------------------------------------------#


# Patient case statements to include descriptive columns
patient_data <- patient_data %>%
  mutate(diagnosis = case_when(
    (`ER STATUS` == 1 | `ER STATUS` == 0) & (`PR STATUS` == 1 | `PR STATUS` == 0) & (`HER2 STATUS` == 1 | `HER2 STATUS` == 0)~ "Triple negative",
    (`ER STATUS` == 2 | `PR STATUS` == 2) & (`HER2 STATUS` == 1 | `HER2 STATUS` == 0) ~ "Luminal A",
    (`ER STATUS` == 2 | `PR STATUS` == 2) & (`HER2 STATUS` == 2 |`HER2 STATUS` == 1 | `HER2 STATUS` == 0) ~ "Luminal B",
    (`ER STATUS` == 1 | `ER STATUS` == 0) & (`PR STATUS` == 1 | `PR STATUS` == 0) & `HER2 STATUS` == 2 ~ "HER2 enriched",
    TRUE ~ "Unknown"
  ))

patient_data <- patient_data %>%
  mutate(Histology_description = case_when(
    HISTOLOGY == 0 ~ "DCIS",
    HISTOLOGY == 1 ~ "Invasive ductal",
    HISTOLOGY == 2 ~ "Lobular",
    HISTOLOGY == 3 ~ "Invasive mucinous",
    HISTOLOGY == 4 ~ "Mixed",
    HISTOLOGY == 5 ~ "Papillary"
  ))

patient_data <- patient_data %>%
  mutate(type = case_when(
    grepl("normal", Full_Filename, ignore.case = TRUE) ~ "Normal",
    grepl("edge", Full_Filename, ignore.case = TRUE) ~ "Edge",
    grepl("tumour", Full_Filename, ignore.case = TRUE) ~ "Tumour",
    TRUE ~ "Undefined"
  ))


patient_data$diagnosis <- factor(patient_data$diagnosis, levels = c("Luminal A","Luminal B","HER2 enriched","Triple negative"))
patient_data$Histology_description <- factor(patient_data$Histology_description, levels = c("DCIS","Invasive ductal","Lobular","Papillary","Mixed"))
patient_data$type <- factor(patient_data$type, levels = c("Normal", "Edge", "Tumour"))

#-------------------------------------------------Filtering--------------------------------------------------------------------#

# Creating Full_filename in cells for parity between formated file names within patient data and the file names in cells. 
cells$Full_Filename <- sub("_\\d+_\\d+_\\d+_\\d+\\.tiff", ".tif", cells$Filename)

# Joining patient and single cell data to allow for easier visualisation
cells <- inner_join(cells, patient_data, by = "Full_Filename")

markers_columns <- c("Nuclear","CD31","CK5","SMA","Ki67","CK8","ANAX1","ECAD","PCNT","CCASP3")

unfiltered_cells <- cells 

# We tried a wide variety of data normalisation techniques. Across all techniques min-max normalisation displayed the best distribution of marker expressions across cells
cells[,markers_columns] <- cells[,markers_columns] %>% mutate(across(everything(), ~ ( . - min(.)) / (max(.) - min(.))* 10))

# Filtering smaller and "lower quality" cells within our dataset. This follows similar methodology to literature
cells <- cells %>% filter(area > 400)
cells <- cells %>% filter(perimeter > 90)

#Renaming columns for IMHC tools
names(cells)[names(cells) == "Filename"] <- "imageNb"
names(cells)[names(cells) == "x"] <- "Pos_X"
names(cells)[names(cells) == "y"] <- "Pos_Y"

names(cells)[names(cells) == "centroid-0"] <- "Pos_X"
names(cells)[names(cells) == "centroid-1"] <- "Pos_Y"

#----------------------------------------------------------------------------------------------------------#
#                                        Classification                                                    #
#----------------------------------------------------------------------------------------------------------#

# Defined cell phenotyping thresolds which were chosen based on characterisation of the original images with their phenotypic classifications by experts

cells <- cells %>%
  mutate(class = case_when(
    (SMA >= 1.2) & (CK8 >= 1)~ "LnB",
    (SMA >= 1.2) & (ECAD >= 0.8)~ "BnL",
    CK8 >= 1 ~ "Luminal",
    ECAD >= 2.5 ~ "Luminal",
    CD31 >= 0.75 ~ "Blood Vessel",
    SMA >= 1 ~ "Basal",
    TRUE ~ "Undefined"
  ))

cells <- cells %>% 
  mutate(proliferating = case_when(
    Ki67 >= 1 ~ "Positive",
    Ki67 < 1 ~ "Negative"
  ))

cells <- cells %>%
  mutate(Type = case_when(
    grepl("normal", imageNb, ignore.case = TRUE) ~ "Normal",
    grepl("edge", imageNb, ignore.case = TRUE) ~ "Edge",
    grepl("tumour", imageNb, ignore.case = TRUE) ~ "Tumour",
    TRUE ~ "Undefined"
  ))

#----------------------------------------------------------------------------------------------------------#
#                                         Cohort Heatmaps                                                  #
#----------------------------------------------------------------------------------------------------------#
unique_filenames <- cells %>%
  select(Full_Filename, Type, diagnosis, Histology_description, GRADE) %>%
  distinct(Full_Filename, .keep_all = TRUE)

heatmap_data <- unique_filenames %>%
  pivot_longer(cols = c(Type, diagnosis, Histology_description, GRADE),
               names_to = "Category",
               values_to = "Value")

ggplot(heatmap_data, aes(x = Full_Filename, y = Category, fill = Value)) +
  geom_tile() +
  scale_fill_manual(values = c("Normal" = "#3274A1", "Edge" = "#FF812C", "Tumour" = "#3A923A",
                               "Luminal A" = "#A9DE00", "Luminal B" = "#60C943", "HER2 enriched" = "#3EA074", "Triple negative" = "#006769",
                               "DCIS" = "#FFBEFF", "Invasive ductal" = "#FF70FF", "Lobular" = "#CB64B8", "Papillary" = "#C964D1", "Mixed" = "#A05EA3",
                               "0" = "#1E90FF","1" = "#4682B4", "2" = "#5F9EA0", "3" = "#B0C4DE")) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 40, vjust = 0.5,size = 7)) +
  labs(x = "Filename", y = "Category", fill = "Classification")

# Plot with ggplot2
ggplot(heatmap_data, aes(x = Full_Filename, y = Category, fill = Value)) +
  geom_tile(color = "white") +  # Adding white border to tiles
  scale_fill_manual(values = color_palette) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 40, vjust = 0.5, size = 7),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Filename", y = "Category", fill = "Classification") +
  geom_rect(data = unique_filenames, 
            aes(xmin = as.numeric(Full_Filename) - 0.5, xmax = as.numeric(Full_Filename) + 0.5, 
                ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1)

#----------------------------------------------------------------------------------------------------------#
#                                        Mean Calculations                                                 #
#----------------------------------------------------------------------------------------------------------#
# For visualisations and further investigation, we calculated the mean values for markers and cell types across our cohort. 


# Transforms columns as factors for ease in plotting. 
cells$Type <- factor(cells$Type, levels = c("Normal", "Edge", "Tumour"))
cells$diagnosis <- factor(cells$diagnosis, levels = c("Luminal A","Luminal B","HER2 enriched","Triple negative"))
cells$Histology_description <- factor(cells$Histology_description, levels = c("DCIS","Invasive ductal","Lobular","Papillary","Mixed"))
cells$GRADE <- factor(cells$GRADE)

cells$class <- factor(cells$class, levels = c("Basal","BnL","Luminal","LnB","Blood Vessel","Undefined"))

total_mean_expr <- cells %>% 
  group_by(imageNb) %>% 
  summarise(across(markers_columns, mean, na.rm = TRUE))

total_mean_expr <- total_mean_expr %>%
  rename_with(~ paste0(.x, "_mean"), markers_columns)

metadata <- cells %>%
  select(imageNb, everything()) %>%
  distinct(imageNb, .keep_all = TRUE)

# Merge mean expressions with metadata
final_df <- metadata %>%
  left_join(total_mean_expr, by = "imageNb")


write.csv(final_df,file = "mean_core_intensity.csv")