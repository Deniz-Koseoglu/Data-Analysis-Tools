#Data Analysis Tools v0.1
#EXAMPLE USE CODE
#Author: Deniz Can Koseoglu
#Date: 01.12.2019

#Source Data Analysis Tools v0.1
source("D:/Data Analysis Tools/Data_Analysis_Tools_v01.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#EXAMPLE 1: PCA with supplementary variables and clustering
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#The data is loaded and analysed using a .CSV file (DATools_PCA_EXAMPLE.csv)
#Create colours for the graph
Color_1 <- rgb(222, 0, 0, maxColorValue=255)
Color_2 <- rgb(231, 151, 45, maxColorValue=255)
Color_3 <- rgb(59, 184, 38, maxColorValue=255)
cols_solid <- c(Color_1, Color_2, Color_3)

#Run PCA
PCA_Results <- Extended_PCA(data = "D:/Data Analysis Tools/DATools_PCA_Example.csv",
                            vars = list(c("Biomarker_1", "Apr", "May", "Jun", "Jul", "Aug", "Apr-May", "SpChla"), 
                                        c("SpChla", "III", "IV", "WD_m", "SpSIC", "SuPIC", "SpPIC", "SpPAR", "SuPAR", "SpSST", "SuSST", "PIP25")),
                            obs = "default", 
                            plot_grphcs = list(c("arrow", "text"), c("point")), 
                            plot_colvar = "contrib",
                            plot_cols = list(cols_solid, c("#00AFBB", "#E7B800", "#FC4E07", brewer.pal(3, "BuGn")), "blue", "grey", "#E7B800"),
                            rank_vars = TRUE, 
                            add_opts = c(TRUE, TRUE, TRUE), 
                            ind_groups = list("HC", "ward.D"),
                            clust_k = c("auto", "auto"),
                            clust_dist = "euclidean",
                            clust_which = c("var", "ind"),
                            export_path = "D:/Data Analysis Tools/Output/Extended_PCA")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#EXAMPLE 2: Agglomerative Hierarchical Clustering (AHC) of time series observations, with accompanying plots
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Load the data
AHC_data <- data.table::fread("D:/Data Analysis Tools/DATools_AHC_EXAMPLE.csv", data.table=FALSE)

#Run AHC
AHC_Results <- Extended_HC(data = AHC_data, 
                           cluster_by="observations", 
                           var_subset = colnames(AHC_data)[-grep("Age", colnames(AHC_data))],
                           var_labels = 1:nrow(AHC_data),
                           dist_measure="euclidean",
                           agglomeration_methods=c("single", "average", "complete", "ward.D", "ward.D2", "mcquitty"),
                           agglomeration_labels=c("Single", "Average", "Complete", "ward.D", "ward.D2", "McQuitty"),
                           preproc_method="zscore",
                           k=3, 
                           cor_method="cophenetic",
                           dataset_label="MD99_2269", 
                           cluster_colours=c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2")),
                           draw_rect=TRUE,
                           coef_method="spearman",
                           export_plots="pdf",
                           export_path="D:/Data Analysis Tools/Output/Extended_HC")

AHC_Overlay <- Overlay_Clusters(data = AHC_data, 
                                cluster_data = AHC_Results[[1]], 
                                cluster_labels = 1:3, 
                                var_subset = colnames(AHC_data)[-grep("Age", colnames(AHC_data))],
                                x_var = "Age", 
                                x_lab = "Age (cal kyr BP)", 
                                colours = cols_solid, 
                                facet_by = "linkage method", 
                                cluster_type = "HC", 
                                k = 3, 
                                export_plots = "pdf",
                                export_path = "D:/Data Analysis Tools/Output/Overlay_Clusters",
                                width = 10,
                                point_size = 9)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#EXAMPLE 3: DTW and alignment of two XRF Ca downcore records to find sensible tie points
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Load data
DTW_list <- lapply(list.files("D:/Data Analysis Tools", full.names = TRUE, pattern = "DTW.*\\.csv$"), data.table::fread, data.table = FALSE)

DTW_Results <- DTW_alignment(query_data = DTW_list[[1]], 
                             ref_data = DTW_list[[2]], 
                             query_title = "GC08", 
                             ref_title = "GS14", 
                             query_subset = colnames(DTW_list[[1]])[-grep("Depth_cm", colnames(DTW_list[[1]]))], 
                             ref_subset = colnames(DTW_list[[2]])[-grep("Depth_cm", colnames(DTW_list[[2]]))],
                             z_norm = list(TRUE, FALSE, FALSE), 
                             reinterp = FALSE,
                             x_var = "Depth_cm",
                             x_label = "Depth (cm)", 
                             y_labels_query = gsub("_ppm", " (ppm)", colnames(DTW_list[[1]])[-grep("Depth_cm", colnames(DTW_list[[1]]))]),
                             y_labels_ref = gsub("_ppm", " (ppm)", colnames(DTW_list[[2]])[-grep("Depth_cm", colnames(DTW_list[[1]]))]),
                             step_pattern = "MVM", 
                             window_type = "sakoechiba", 
                             window_size = 500, 
                             MVM_elasticity = 50, 
                             RJ_step_settings = c(4, "d"), 
                             open_begin = FALSE, 
                             open_end = FALSE, 
                             x_align = TRUE,
                             y_offset = c(1.5, "new_breaks"),
                             match_min = 20, 
                             sample_freq = 1, 
                             match_vis = 0.8, 
                             grang_order = 3,
                             x_rounding_order = -1, 
                             y_rounding_order = -1,
                             dtw_distance = "Euclidean",
                             export_plots = "pdf", 
                             export_path = "D:/Data Analysis Tools/Output/DTW_alignment")
