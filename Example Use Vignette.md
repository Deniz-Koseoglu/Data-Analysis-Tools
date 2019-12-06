# Using Data Analysis Tools
This vignette is designed to showcase the functionality of Data Analysis Tools v0.1 (December 2019) and assumes you have [downloaded and installed](https://github.com/Deniz-Koseoglu/Data-Analysis-Tools#getting-started) the necessary [code](https://github.com/Deniz-Koseoglu/Data-Analysis-Tools/blob/master/Data_Analysis_Tools_v01.R) and [example data](https://github.com/Deniz-Koseoglu/Data-Analysis-Tools/blob/master/DATools_EXAMPLE_DATA.zip). A basic understanding of R syntax (v3.5.0) is also advantageous. The vignette assumes that both the source code and example data were unpacked to **D:/Data Analysis Tools**. All examples shown were are also available in the [accompanying .R file](https://github.com/Deniz-Koseoglu/Data-Analysis-Tools/blob/master/Data_Analysis_Tools_v01_EXAMPLES.R).

### Example 1: PCA with clustering
The `Extended_PCA` function is used here to:
1. Assess the clustering tendency of data using a visual heatmap-based technique and the Hopkins' Statistic.
2. Reduce the dimensionality of a multivariate dataset to achieve a 2-dimensional representation.
3. Cluster the resulting observation and variable PC coordinates to identify any significant groupings.
4. Visualize the results for all of the above using `ggplot2`.

Run the following code:
```r
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
```

Here are some examples of graphical .PDF output, saved to **D:/Data Analysis Tools/Output/Extended_PCA** in this case:
The eigenvalue bar plot reveals that the first two PCs explain <74% of the variance attributable to primary variables/observations.
![Image1](https://i.ibb.co/1LwmzxF/Easy-PCA-Plots-2019-12-03-05hr-46min-54sec-Page-01.png)
The correlation circle (for primary and supplementary variables) is shown below, along with primary observation coordinates along the first two PCs. Both primary variables and observations are separated into two clusters by Ward's-linkage AHC, where the optimal cluster number was determined automatically using Cluster Validity Indices.
![Image2](https://i.ibb.co/TR2nrQ8/Easy-PCA-Plots-2019-12-03-05hr-46min-54sec-Page-04.png)
The variable plot shows that one distinct cluster encompassing variables **Jun, Jul, Aug**, which contributes primarily to PC2. Along PC1, we see variables **Biomarker_1, April, May, April-May, April-June** strongly associated with each other and described primarily by PC1. The supplementary variables are shown by dashed blue lines, exhibit much smaller magnitude relative to primary variables, and are thus not very well-represented on the PC map.
![Image3](https://i.ibb.co/ZScmQKX/Easy-PCA-Plots-2019-12-03-05hr-46min-54sec-Page-05.png)
The plots of observations/individuals shows two clusters that are mainly separated along PC1 (x-axis) but are independent of PC2 (both exhibit positive and negative coordinates along the y-axis). The wide, near-uniform distribution of PC coordinates, together with overlapping ellipses around centroids suggests that the observations are probably not readily-clusterable in this case - at least not along the first two PCs! Using a density-based clustering approach (e.g. DBSCAN) in thise case would likely yield poor results.
Contributions of both variables and observations to the first two PCs may also be visualised using bar charts, as below:
![Image4](https://i.ibb.co/RCXp3Yv/Easy-PCA-Plots-2019-12-03-05hr-46min-54sec-Page-07.png)
![Image5](https://i.ibb.co/fCJ0P3Z/Easy-PCA-Plots-2019-12-03-05hr-46min-54sec-Page-10.png)
Variable eigenvalues, variable and observation coordinates and squared cosines across the first 5 PCs (since `pc_keep = 5` in this instance), cluster memberships and other results are also available as raw data in the two .CSV output files accompanying .PDF plots.
<br></br>
### Example 2: AHC time series clustering and visualisation of results
This example uses the `Extended_HC` and `Overlay_Clusters` functions to cluster a 4-dimensional time series as follows:
```r
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
                                facet_by = "variables", 
                                cluster_type = "HC", 
                                k = 3, 
                                export_plots = "pdf",
                                export_path = "D:/Data Analysis Tools/Output/Overlay_Clusters",
                                width = 5,
                                point_size = 9)
```
Example output below shows the AHC dendrogram output using `ggplot2` and Ward's linkage method:
![Image6](https://i.ibb.co/c6672yF/HC-Clustering-by-Observations-2019-12-03-05hr-47min-30sec-Page-21.png)
The label size appears to be too large to display all observations clearly! This can be adjusted, where necessary, via the `gg_lab_size` argument. Alternatively, `var_labels` may be set to `NULL` or `NA` to remove them entirely.
A cophenetic Spearman correlation matrix between all utilized linkage methods is shown below. It is apparent that single linkage is most dissimilar from the others, while average and complete linkage are also dissimilar to the **ward.D** method.
![Image8](https://i.ibb.co/JrMwskq/HC-Clustering-by-Observations-2019-12-03-05hr-47min-30sec-Page-23.png)
The accompanying .CSV output files also include the results of Cluster Validity Indices (CVI) to determine a suitable number of clusters, as well as cluster memberships using each linkage method.
Finally, the output of `Overlay_Clusters` using the results of `Extended_HC` yields the following plots of cluster memberships (the example shown uses Ward's linkage method) overlaid onto each of the four variables included in the dataset:
![Image9](https://i.ibb.co/g9KjPdp/Cluster-Overlay-by-Variable-2019-12-05-19hr-46min-06sec-Page-4.png)
In this instance, Ward's linkage AHC clustering effectively identified the increase in Variable_4, in particular. 
<br></br>
### Example 3: DTW alignment of time series from two downcore records
Differences in sedimentation rates contribute to the differences observed between depth-age profiles of downcore sedimentary records. Radiocarbon dating of carbonate deposited within ancient sediment is costly, time-consuming, and only possible where sufficient preserved material is available. Thus, even downcore records considered well-dated often have multi-centennial age resolution at maximum. However, where at least two cores from the same location are available, their age models may be merged by finding tie-points between the same measured quantity(ies) across a range of core depths (e.g. calcium and other trace metals); such tie-points often take the form of significantly similar shape-based features, such as maxima and minima. Thus, the depth (and, therefore, age) of one core may be inferred from another based on these similarities. Here, we will use `DTW_alignment` to identify tie-points from two XRF-measured calcium depth profiles obtained from equi-located sediment cores.
