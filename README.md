# Data Analysis Tools
A set of statistical workflows implemented in [RStudio](https://www.rstudio.com/products/rstudio/download/) (R v3.5.0) and mostly consisting of wrappers for dimensionality reduction (e.g. PCA), unsupervised classification (e.g. clustering), and other data analysis functions.

Please send any suggestions for bug fixes and any other queries to deniz.koseoglu@plymouth.ac.uk.

**NOTE**: Textual (.CSV) and visual (.PDF/.PNG) output produced by Data Analysis Tools v0.1, and toolkit functionality, are showcased in the accompanying [vignette]().

Full functionality of the tools was confirmed only using Windows 10. Data Analysis Tools are provided as-is under MIT Licence terms (see LICENCE.md in the repository root).

# Functionality
Data Analysis Tools v0.1 include the following functions:
1. `Extended_PCA` determines clustering tendency of input data, carries out PCA with or without prior data scaling and inclusion of supplementary variables/observations. Additionally, both variables and/or observation PCA coordinates may be clustered using a selection of different, adjustable methodologies, with or without automatic determination of an optimal cluster number.
2. `Extended_HC` visualises clustering tendency, as well as performs and visualises Agglomerative Hierarchical Clustering (AHC) of input data using a selection of 6 linkage methods. 
3. `Overlay_Clusters` plots cluster memberships overlaid onto a scatter plot. This is useful after clustering timeseries by observations/rows (rather than variables/columns). A specialised method for the results of `Extended_HC` exists for visualising data, but any factor of cluster memberships may be passed to the function.
4. `DTW_alignment` performs Dynamic Time Warping (DTW) on any number of paired time series, which may be of different length. The objective is to align the time series based on their similarity. The time series may be of different length and timestep, of which the latter may also be irregularly spaced. A DTW distance matrix between the time series is produced, allowing out-of-phase time series to be optimally aligned, provided they share identifiable similarities in shape. Thus, the utility of DTW is apparent when determining tie-points between paleoclimatological time series, for example.

# Dependencies
R (≥3.5.0), [biotools](https://cran.r-project.org/web/packages/biotools/index.html), [asbio](https://cran.r-project.org/web/packages/asbio/index.html), [readr](https://cran.r-project.org/web/packages/readr/index.html), [lattice](https://cran.r-project.org/web/packages/lattice/index.html), [corrplot](https://cran.r-project.org/web/packages/corrplot/index.html), [fmsb](https://cran.r-project.org/web/packages/fmsb/index.html), [latticeExtra](https://cran.r-project.org/web/packages/latticeExtra/index.html), [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html), [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html), [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html), [car](https://cran.r-project.org/web/packages/car/index.html), [psych](https://cran.r-project.org/web/packages/psych/index.html), [data.table](https://cran.r-project.org/web/packages/data.table/index.html), [gridGraphics](https://cran.r-project.org/web/packages/gridGraphics/index.html), [ggcorrplot](https://cran.r-project.org/web/packages/ggcorrplot/index.html), [corrgram](https://cran.r-project.org/web/packages/corrgram/index.html), [Hmisc](https://cran.r-project.org/web/packages/Hmisc/index.html), [xlsx](https://cran.r-project.org/web/packages/xlsx/index.html), [factoextra](https://cran.r-project.org/web/packages/factoextra/index.html), [NbClust](https://cran.r-project.org/web/packages/NbClust/index.html), [dendextend](https://cran.r-project.org/web/packages/dendextend/index.html), [fpc](https://cran.r-project.org/web/packages/fpc/index.html), [dtw](https://cran.r-project.org/web/packages/dtw/index.html), [dtwclust](https://cran.r-project.org/web/packages/dtwclust/index.html), [TSclust](https://cran.r-project.org/web/packages/TSclust/index.html), [rlist](https://cran.r-project.org/web/packages/rlist/index.html), [pryr](https://cran.r-project.org/web/packages/pryr/index.html), [GGally](https://cran.r-project.org/web/packages/GGally/index.html), [ggpmisc](https://cran.r-project.org/web/packages/ggpmisc/index.html), [lmtest](https://cran.r-project.org/web/packages/lmtest/index.html), [cpm](https://cran.r-project.org/web/packages/cpm/index.html), [ecp](https://cran.r-project.org/web/packages/ecp/index.html), [changepoint](https://cran.r-project.org/web/packages/changepoint/index.html), [magrittr](https://cran.r-project.org/web/packages/magrittr/index.html), [scales](https://cran.r-project.org/web/packages/scales/index.html).

# Getting Started
1. Make sure you have [R](https://cran.r-project.org/mirrors.html) and [RStudio](https://www.rstudio.com/products/rstudio/download/) installed.
2. Download the [Data Analysis Tools]() source code and place it in a directory of your choice. Examples in this manual and the [vignette]() assume the .R file is located at **D:/Data Analysis Tools**.
3. If you wish, you can also download various [example data files]() used in the [vignette]().
4. Create a new script in RStudio using **File -> New File -> R Script**.
5. Source Data Analysis Tools using the command: 
```r 
source("D:/Data Analysis Tools/Data_Analysis_Tools_v01.R"). 
```
This file may also be opened in RStudio directly to alter the source code, if required.

# Functions
### The `Extended_PCA` function
#### Usage
```r
Extended_PCA <- function(data, vars, obs="default", plot_grphcs=list(c("arrow", "text"), c("point")), plot_colvar="red",
                         plot_cols, rank_vars=TRUE, add_opts=c(TRUE, FALSE, TRUE), ind_groups=c("HC", "complete"), 
                         clust_k="auto", clust_dist="euclidean", clust_which="none", ellipses=TRUE, pc_keep = 5,
                         export_path=getwd(), export_res="all", height=10, width=10, point_size=12)
```

#### Arguments
| Argument | Description |
| ------------- |-------------|
|**arg**|desc.|

#### Details
Please refer to the [Data Analysis Tools vignette]() for example usage and output of this and other functions.

#### Values
<br></br>
### The `Extended_HC` function
#### Usage
```r
Extended_HC <- function(data, cluster_by="observations", var_subset, var_labels=var_subset, dist_measure="euclidean",
                        agglomeration_methods=c("single", "average", "complete", "ward.D", "ward.D2", "mcquitty"),
                        agglomeration_labels=c("Single", "Average", "Complete", "ward.D", "ward.D2", "McQuitty"),
                        preproc_method="zscore", cvi_range=c(2,10), cluster_legend=TRUE, cluster_labels=NULL, k=3,
                        cor_method="cophenetic", horizont=TRUE, dataset_label="Sample", cluster_colours=c(brewer.pal(9, "Set1"),
                        brewer.pal(8, "Set2")), gg_lab_size=0.9, base_lab_size=0.8, gg_lwd=0.8, base_lwd=1.2, draw_rect=TRUE,
                        coef_method="spearman", FM_test=NULL, FM_lim=c(0.7,1), export_plots="pdf", width=5, height=5,
                        point_size=10, dpi=500, export_results=TRUE, export_path=getwd(), dend_mar=c(4,4))
```

#### Arguments
| Argument | Description |
| ------------- |-------------|
|**arg**|desc.|

#### Details
Please refer to the [Data Analysis Tools vignette]() for example usage and output of this and other functions.

#### Values
<br></br>
### The `Overlay_Clusters` function
#### Usage
```r
Overlay_Clusters <- function(data, cluster_data, cluster_labels, var_subset, var_labels=var_subset, x_var, x_lab="X label",
colours, facet_by="linkage method", cluster_type="arbitrary", k, export_plots="pdf", export_results=TRUE, export_path=getwd(), height=5, width=5, point_size=10, dpi=500)
```

#### Arguments
| Argument | Description |
| ------------- |-------------|
|**arg**|desc.|

#### Details
Please refer to the [Data Analysis Tools vignette]() for example usage and output of this and other functions.

#### Values
<br></br>
### The `DTW_alignment` function
#### Usage
```r
DTW_alignment <- function(query_data, ref_data, query_title="Query", ref_title="Reference", query_subset=NULL, ref_subset=NULL,
z_norm=list(TRUE, FALSE, FALSE), reinterp=FALSE, x_var, step_pattern="symmetric2", window_type="sakoechiba", window_size=20,
MVM_elasticity=50, RJ_step_settings=c(4, "d"), open_begin=FALSE, open_end=FALSE, x_align=TRUE, y_offset=NULL, x_label,
y_labels_query=waiver(), y_labels_ref=waiver(), colours=brewer.pal(9, "Set1"), match_subset=NULL, match_min=NULL, 
sample_freq=1, match_vis=0.8, grang_order=3, x_rounding_order=-1, y_rounding_order=-1, export_plots="pdf", width=10, height=5,
point_size=10, dpi=500, export_results=TRUE, export_path=getwd(), dtw_distance="Euclidean")
```

#### Arguments
| Argument | Description |
| ------------- |-------------|
|**arg**|desc.|

#### Details
Please refer to the [Data Analysis Tools vignette]() for example usage and output of this and other functions.

#### Values
