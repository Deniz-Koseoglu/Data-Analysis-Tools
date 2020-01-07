# Data Analysis Tools
A set of statistical workflows implemented in [RStudio](https://www.rstudio.com/products/rstudio/download/) ([R](https://cran.r-project.org/mirrors.html) v3.5.0) and mostly consisting of wrappers for dimensionality reduction (e.g. PCA), unsupervised classification (e.g. clustering), and other data analysis functions.

Please send any suggestions for bug fixes and any other queries to deniz.koseoglu@plymouth.ac.uk.

**NOTE**: Textual (.CSV) and visual (.PDF/.PNG) output produced by Data Analysis Tools v0.1, and toolkit functionality, are showcased in the accompanying [vignette](https://github.com/Deniz-Koseoglu/Data-Analysis-Tools/blob/master/Example%20Use%20Vignette.md).

Full functionality of the tools was confirmed only using Windows 10. Data Analysis Tools are provided as-is under MIT Licence terms (see LICENCE.md in the repository root).

# Functionality
Data Analysis Tools v0.1 include the following functions:
1. `Extended_PCA` determines clustering tendency of input data, carries out PCA with or without prior data scaling and inclusion of supplementary variables/observations. Additionally, both variables and/or observation PCA coordinates may be clustered using a selection of different, adjustable methodologies, with or without automatic determination of an optimal cluster number.
2. `Extended_HC` visualises clustering tendency, as well as performs and visualises Agglomerative Hierarchical Clustering (AHC) of input data using a selection of 6 linkage methods. 
3. `Overlay_Clusters` plots cluster memberships overlaid onto a scatter plot. This is useful after clustering timeseries by observations/rows (rather than variables/columns). A specialised method for the results of `Extended_HC` exists for visualising data, but any factor of cluster memberships may be passed to the function.
4. `DTW_alignment` performs Dynamic Time Warping (DTW) on any number of paired time series, which may be of different length. The objective is to align the time series based on their similarity. The time series may be of different length and timestep, of which the latter may also be irregularly spaced. A DTW distance matrix between the time series is produced, allowing out-of-phase time series to be optimally aligned, provided they share identifiable similarities in shape. Thus, the utility of DTW is apparent when determining tie-points between paleoclimatological time series, for example.

# Dependencies
R (≥3.5.0), [biotools](https://cran.r-project.org/web/packages/biotools/index.html), [asbio](https://cran.r-project.org/web/packages/asbio/index.html), [readr](https://cran.r-project.org/web/packages/readr/index.html), [lattice](https://cran.r-project.org/web/packages/lattice/index.html), [corrplot](https://cran.r-project.org/web/packages/corrplot/index.html), [fmsb](https://cran.r-project.org/web/packages/fmsb/index.html), [latticeExtra](https://cran.r-project.org/web/packages/latticeExtra/index.html), [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html), [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html), [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html), [car](https://cran.r-project.org/web/packages/car/index.html), [psych](https://cran.r-project.org/web/packages/psych/index.html), [data.table](https://cran.r-project.org/web/packages/data.table/index.html), [gridGraphics](https://cran.r-project.org/web/packages/gridGraphics/index.html), [ggcorrplot](https://cran.r-project.org/web/packages/ggcorrplot/index.html), [corrgram](https://cran.r-project.org/web/packages/corrgram/index.html), [Hmisc](https://cran.r-project.org/web/packages/Hmisc/index.html), [xlsx](https://cran.r-project.org/web/packages/xlsx/index.html), [factoextra](https://cran.r-project.org/web/packages/factoextra/index.html), [NbClust](https://cran.r-project.org/web/packages/NbClust/index.html), [dendextend](https://cran.r-project.org/web/packages/dendextend/index.html), [fpc](https://cran.r-project.org/web/packages/fpc/index.html), [dtw](https://cran.r-project.org/web/packages/dtw/index.html), [dtwclust](https://cran.r-project.org/web/packages/dtwclust/index.html), [TSclust](https://cran.r-project.org/web/packages/TSclust/index.html), [rlist](https://cran.r-project.org/web/packages/rlist/index.html), [pryr](https://cran.r-project.org/web/packages/pryr/index.html), [GGally](https://cran.r-project.org/web/packages/GGally/index.html), [ggpmisc](https://cran.r-project.org/web/packages/ggpmisc/index.html), [lmtest](https://cran.r-project.org/web/packages/lmtest/index.html), [cpm](https://cran.r-project.org/web/packages/cpm/index.html), [ecp](https://cran.r-project.org/web/packages/ecp/index.html), [changepoint](https://cran.r-project.org/web/packages/changepoint/index.html), [magrittr](https://cran.r-project.org/web/packages/magrittr/index.html), [scales](https://cran.r-project.org/web/packages/scales/index.html), [BBmisc](https://cran.r-project.org/web/packages/BBmisc/index.html).

# Getting Started
1. Make sure you have [R](https://cran.r-project.org/mirrors.html) and [RStudio](https://www.rstudio.com/products/rstudio/download/) installed.
2. Download the [Data Analysis Tools](https://github.com/Deniz-Koseoglu/Data-Analysis-Tools/blob/master/Data_Analysis_Tools_v01.R) source code and place it in a directory of your choice. Examples in this manual and the [vignette](https://github.com/Deniz-Koseoglu/Data-Analysis-Tools/blob/master/Example%20Use%20Vignette.md) assume the .R file is located at **D:/Data Analysis Tools**.
3. If you wish, you can also download various [example data files](https://github.com/Deniz-Koseoglu/Data-Analysis-Tools/blob/master/DATools_EXAMPLE_DATA.zip) used in the [vignette](https://github.com/Deniz-Koseoglu/Data-Analysis-Tools/blob/master/Example%20Use%20Vignette.md).
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
|**data**|Character filepath to a .CSV file, or the name of an R `data.frame` or `matrix` containing input data.|
|**vars**|A character vector of `data` column names to use as *primary variables* for PCA. If a character `list` of length 2 is provided, the **first** and **second** elements denote *primary* and *secondary/supplementary* variables, respectively.|
|**obs**|The `"default"` value uses all observations/rows for analysis. If a numeric vector of row numbers is provided, these are treated as *primary observations*; *secondary/supplementary observations* can be specified as the **second element** of a numeric list of length 2.|
|**plot_grphcs**|A character vector list of length 2, denoting the plotting geometry of variables and observations, respectively. Possible character values include: `"points"` (data markers), `"arrows"` (arrows from the origin), `"text"` (data labels), or any combination of these. By default, arrows and labels are plotted for variables, while only points are shown for observations.|
|**plot_colvar**|How to colour variables in PCA plots? Values of `"contrib"` and `"cos2"` colour variable PC loadings based on their contribution and squared cosine values, respectively. A single colour for variables can also be specified, e.g. `"red"` (default), but **note** that the PCA biplot containing both variables and observations will not be produced in this case.|
|**plot_cols**|A value or vector list of length 2 specifying variable and observation colours to use, respectively.|
|**rank_vars**|`TRUE/FALSE` logical. Should the correlation circle for variables be derived from their ranks? Essentially, when `TRUE`, Spearman rather than the Pearson correlation is performed on `vars`.|
|**add_opts**|A `TRUE/FALSE` logical vector of length 3. The **first element** determines whether `data` is scaled (converted to z-scores) prior to PCA. The **second** and **third** elements are aesthetic options specifying whether observation cluster centroids should be shown in plots, and whether label text should be repelled, respectively.|
|**ind_groups**|A list specifying the type and settings of clustering method applied to variable and/or observation PC coordinates. The **first element** denotes the clustering method and is one of: `"HC"` (Agglomerative Hierarchical Clustering, AHC; default), `"dbscan"` (DBSCAN algorithm), `"fcm"` (fuzzy c-means). A value of `FALSE` foregoes clustering. The **second element** is either the linkage method when AHC is used (one of: `"single"`, `"average"`, `"complete"` (default), `"ward.D"`, `ward.D2`, `"mcquitty"`), or a numberic value for epsilon neighbourhood size (EPS) when using DBSCAN. The **third element** is only used for DBSCAN and is an integer specifying the minimum number of nearest neighbours each identified cluster should have.|
|**clust_which**|Any of: `"var"` and/or `"ind"`. Specifies whether variables and/or observations should be clustered, respectively.|
|**clust_k**|Either a numeric value(s) specifying the number of clusters to separate variables and/or observations into (as specified in `clust_which`), or the value `"auto"`, which attempts to determine an optimal cluster number(s) automatically. If a vector of length 2 is provided, the **first** and **second** element are used for clustering variables and observations, respectively. Alternatively, if a single value is provided, it is either propagated to both variables and observations.|
|**clust_dist**|A character value of the distance measure to use. Defaults to `"euclidean"`. See `?fviz_dist` (`method` argument) for available values.|
|**ellipses**|`TRUE/FALSE` logical. Should centroid ellipses be drawn around clustered observations?|
|**pc_keep**|A numeric value specifying the number of the PCs to keep (5 by default).|
|**export_path**|Character path to the directory where results (including plots) should be exported. Defaults to the working directory.|
|**height**|The height of exported .PDF graphics (10 by default).|
|**width**|The width of exported .PDF graphics (10 by default).|
|**point_size**|A numeric value of point size to use in exported .PDF graphics (12 by default).|

#### Details
Please refer to the [Data Analysis Tools vignette](https://github.com/Deniz-Koseoglu/Data-Analysis-Tools/blob/master/Example%20Use%20Vignette.md) for example usage and output of this and other functions.

#### Values
Returns a list of length 2, where the **first element** is an object of class `PCA` and contains comprehensive results of the analysis. The **second element** includes all plot objects produced by the function via `ggplot2`.
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
|**data**|The name of an R `data.frame` or `matrix` containing input data.|
|**cluster_by**|One of `"variables"` or `"observations"`. Denotes whether `data` should be clustered by columns or rows, respectively.|
|**var_subset**|A character vector of column names from `data` to use.|
|**var_labels**|An optional custom vector of labels for variables specified in `var_subset`.|
|**dist_measure**|The distance measure to use. Defaults to `"euclidean"`. Available options are outlined in `?proxy::dist` and `?tsclust`.|
|**agglomeration_methods**|Character value/vector of linkage methods to use. Any of: `"single"`, `"average"`, `"complete"`, `"ward.D"`, `"ward.D2"`, `"mcquitty"`. All linkage methods are used by default.|
|**agglomeration_labels**|Character vector/value of plot labels for linkage methods specified in `agglomeration_methods`.|
|**preproc_method**|Character. The type of pre-processing applied to `data`. Currently supported values are either `NULL` or `"zscore"`.|
|**cvi_range**|A numeric vector of length 2, where the **first** and **second** elements represent the minimum and maximum number of clusters to derive Cluster Validity Indices for.|
|**cluster_legend**|`TRUE/FALSE` logical. Should a cluster legend be displayed in plots?|
|**cluster_labels**|An optional vector of labels for `cluster_legend`.|
|**k**|The number of clusters to split dendrograms into (3 by default).|
|**cor_method**|A character value giving the method to use when deriving a correlation matrix for comparison of different AHC linkage types (as determined by `agglomeration_methods`). One of: `"cophenetic"` (default), `"baker"`, `"common_nodes"`, or `"FM_index"`.|
|**coef_method**|The type of correlation coefficient to derive for `cor_method`. One of: `"spearman"` (default), `"pearson"`, or `"kendall"`.|
|**horizont**|`TRUE/FALSE` logical. Should dendrograms be oriented horizontally?|
|**dataset_label**|A character value containing an optional `data` label/title.|
|**cluster_colours**|A numeric or character vector of colours to use for clusters.|
|**gg_lab_size**|Numeric `ggplot2` label text size (0.9 by default).|
|**base_lab**|Numeric R base graphics label text size (0.8 by default).|
|**gg_lwd**|Numeric `ggplot2` line width (0.8 by default).|
|**base_lwd**|Numeric R base graphics line width (1.2 by default).|
|**draw_rect**|A `TRUE/FALSE` logical. Should dashed rectangles be drawn around identified clusters within dendrograms?|
|**FM_test**|`NULL` by default. Otherwise, a column name containing cluster/class memberships for observations **previously** specified by the user. These are used to carry out Fowlkes-Mallows testing of said membeships against those determined via specified `agglomeration_methods`. Results showing the performance of each linkage method are subsequently visualised in a graph. Testing is only carried out when `cluster_by = "observations"`.|
|**FM_lim**|A numeric vector of length 2 and values between 0 and 1. Determines the axis limits for the Fowlkes-Mallows testing visualisation.|
|**export_plots**|One of: `"none"`, `"png"`, or `"pdf"` (default). Specifies whether plots are exported, and the filetype.|
|**width**|The width of exported .PDF/.PNG graphics (10 by default).|
|**height**|The height of exported .PDF/.PNG graphics (10 by default).|
|**point_size**|The point size of exported .PDF/.PNG graphics (10 by default).|
|**dpi**|The resolution of exported .PNG plots (when `export_plots = "png"`), in Dots Per Inch (500 by default).|
|**export_results**|`TRUE/FALSE` logical. Should results be exported as a .CSV file?|
|**export_path**|Character path to the directory where results (including plots) should be exported. Defaults to the working directory.|
|**dend_mar**|A numeric vector of length 2 providing the bottom (**first element**) and right (**second element**) margins for dendrograms drawn using base R graphics (rather than `ggplot2`).|

#### Details
Please refer to the [Data Analysis Tools vignette](https://github.com/Deniz-Koseoglu/Data-Analysis-Tools/blob/master/Example%20Use%20Vignette.md) for example usage and output of this and other functions.

#### Values
Analogously to `Extended_PCA`, the function returns a list of length 2, where the **first element** contains clustering results for each of chosen `agglomeration_methods`, and the **second element** contains graphical `ggplot2` output.
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
|**data**|The name of an R `data.frame` or `matrix` containing input data.|
|**cluster_data**|Clustering information for `data`. Either a `HierarchicalTSClusters` object (e.g. from `Extended_HC` or `tsclust` output), or a factor vector of cluster memberships.|
|**cluster_labels**|An optional character vector of cluster labels.|
|**var_subset**|Character vector of names for `data` columns to be plotted.|
|**var_labels**|An optional character vector of custom `var_subset` labels.|
|**x_var**|The column name within `data` containing the variable to be plotted on the x-axis.|
|**x_lab**|Optional custom character label for `x_var`.|
|**colours**|Character or numeric vector of colours to use for clusters in plots.|
|**facet_by**|One of: `"linkage_method"` (default), `"variables"`, or `"none"`. Should output plots be faceted (by either variable or linkage method, when `cluster_type = "HC"`)?|
|**cluster_type**|One of: `"HC"` or `"arbitrary"` (default), depending on whether the class of `cluster_data` is `HierarchicalTSClusters` (e.g. the output of `Extended_HC`) or a simple factor of cluster memberships, respectively.|
|**k**|The number of clusters to cut dendrograms into. **NOTE** that this only applies when `cluster_type = "HC"`.|
|**export_plots**|One of: `"none"`, `"png"`, or `"pdf"` (default). Specifies whether plots are exported, and the filetype.|
|**export_results**|`TRUE/FALSE` logical. Should results be exported as a .CSV file?|
|**export_path**|Character path to the directory where results (including plots) should be exported. Defaults to the working directory.|
|**height**|The height of exported .PDF graphics (5 by default).|
|**width**|The width of exported .PDF graphics (5 by default).|
|**point_size**|The point size of exported .PDF graphics (10 by default).|
|**dpi**|The resolution of exported .PNG plots (when `export_plots = "png"`), in Dots Per Inch (500 by default).|

#### Details
Please refer to the [Data Analysis Tools vignette](https://github.com/Deniz-Koseoglu/Data-Analysis-Tools/blob/master/Example%20Use%20Vignette.md) for example usage and output of this and other functions.

#### Values
Currently returns a list of `ggplot2` plots created.
<br></br>
### The `DTW_alignment` function
#### Usage
```r
DTW_alignment <- function(query_data, ref_data, query_title="Query", ref_title="Reference", query_subset=NULL, ref_subset=NULL,
z_norm=list(TRUE, FALSE, FALSE), reinterp=FALSE, x_var, step_pattern="symmetric2", window_type="sakoechiba", window_size=20,
MVM_elasticity=50, RJ_step_settings=c(4, "d", FALSE), open_begin=FALSE, open_end=FALSE, x_align=TRUE, y_offset=NULL, x_label=x_var,
y_labels_query=waiver(), y_labels_ref=waiver(), colours=brewer.pal(9, "Set1"), match_subset=NULL, match_min=NULL, 
sample_freq=1, match_vis=0.8, grang_order=3, x_rounding_order=-1, y_rounding_order=-1, export_plots="pdf", width=10, height=5,
point_size=10, dpi=500, export_results=TRUE, export_path=getwd(), dtw_distance="Euclidean")
```

#### Arguments
| Argument | Description |
| ------------- |-------------|
|**query_data**|A `data.frame` or `matrix` object to use as the query dataset (to be matched via DTW to `ref_data`).|
|**ref_data**|A `data.frame` or `matrix` object to use as the reference dataset (to which `query_data` is matched).|
|**query_title**|An optional title for `query_data`.|
|**ref_title**|An optional title for `ref_data`.|
|**query_subset**|A character vector of column names from `query_data` to use for analysis. Paired element-wise to `ref_subset`.|
|**ref_subset**|A character vector of column names from `ref_data` to use for analysis. Paired element-wise to `query_subset`.|
|**z_norm**|A list of length 3. The **first element** is a `TRUE/FALSE` logical; when `TRUE`, variable data are converted to z-scores prior to DTW alignment. The **second** and **third** elements are either `FALSE` or numeric vectors of length 2, each denoting the minimum and maximum values to which `query_data` and `ref_data` should be normalized within plots, respectively. For example, if `z_norm[[2]] == c(0,1)`, `query_data` will be normalized to between 0 and 1. **NOTE**: When such normalisation is carried out (i.e. `z_norm[[2]]` and `z_norm[[3]]` are *not* `FALSE`), this overrides `y_offset`!|
|**reinterp**|When `TRUE`, linearly reinterpolates `query_data` to the length of `ref_data` prior to analysis. `FALSE` by default. Reinterpolation should not cause a significant difference to the analysis; for example, see [this article](http://www.cs.ucr.edu/~eamonn/DTW_myths.pdf).|
|**x_var**|Character column name containing variable data to plot on the x-axis (e.g. core depth, age).|
|**x_label**|Character. Name for `x_var` (column name is used by default).|
|**step_pattern**|Character value representing the step pattern to use. See the `step.pattern` argument for `?dtw`. Currently supported values: `symmetric1`, `symmetric2` (default), `asymmetric`, `RJ` (Rabiner-Juang; see `?rabinerJuangStepPattern`), or `MVM` (Minimum Variance Matching algorithm; see `?mvmStepPattern`).|
|**window_type**|A character value denoting the type of DTW windowing function. One of: `"none"`, `"itakura"`, `"sakoechiba"` (default), or `"slantedband"`. See also `?dtw`.|
|**window_size**|The numeric value denoting the maximum permitted size of the warping window.|
|**MVM_elasticity**|Elasticity (numeric value) when `step_pattern = "MVM"`.|
|**RJ_step_settings**|A vector of length 3 providing slope weighting, path specification, and smoothing when `step_pattern = "RJ"`. See `?rabinerJuangStepPattern` for details.|
|**open_end**|When `TRUE`, performs an open-ended alignment, whereby matching of `query_data` in its entirety to a **leading** part of `ref_data` is carried out.|
|**open_begin**|When `TRUE` together with `open_end`, matches all of `query_data` to any continuous part of `ref_data`, without a pre-defined beginning or end. Useful for sub-sequence finding.|
|**x_align**|`TRUE/FALSE` logical. Should `query_data` timeseries specified in `query_subset` be visually matched in length to those provided in `ref_data` and `ref_subset`? Works well for time series with significantly different x-axis values.|
|**y_offset**|A numeric value, or vector of length 2. The **first_element** determines the extent to which `query_data` is offset along the y-axis from `ref_data` in output plots. This is defined through multiplication of `y_offset[1]` with the arithmetic mean of `query_data`. Defaults to `NULL` (no offset). Optionally, `"new_breaks"` can be specified as the **second element**; this will create a second set of axis breaks (and a secondary y-axis) for the newly offset `query_data`. Setting a `y_offset` works well for time series with significantly overlapping y-axis values.|
|**sample_freq**|The frequency with which to draw DTW lines between `query_data` and `ref_data` in output graphics. Defaults to `1` (a line is drawn for every available sample).|
|**match_vis**|A number from 0 to 1 specifying the visibility of DTW lines connecting timeseries (0.8, i.e. 80% visibility, by default).|
|**grang_order**|Numeric value of the lag order to use when testing for Granger causality between `query_data` and `ref_data`. See also `?grangertest`.|
|**x_rounding_order**|Rounding order to use when deriving x-axis breaks for ggplot2 graphics. Defaults to `-1` (rounded to the nearest 10; see also `?round`).|
|**y_rounding_order**|Rounding order to use when deriving y-axis breaks for ggplot2 graphics. Defaults to `-1` (rounded to the nearest 10; see also `?round`).|
|**export_plots**|One of: `"none"`, `"png"`, or `"pdf"` (default). Specifies whether plots are exported, and the filetype.|
|**width**|The width of exported .PDF/.PNG plots (10 by default).|
|**height**|The height of exported .PDF/.PNG plots (5 by default).|
|**point_size**|The point size of exported .PDF/.PNG plots (10 by default).|
|**dpi**|The resolution of exported .PNG plots (when `export_plots = "png"`), in Dots Per Inch (500 by default).|
|**export_results**|`TRUE/FALSE` logical. Should results be exported as a .CSV file?|
|**export_path**|Character path to the directory where results (including plots) should be exported. Defaults to the working directory.|
|**dtw_distance**|**Currently not utilized**. Character value denoting the pointwise (local) distance to use for multivariate timeseries.|

#### Details
Please refer to the [Data Analysis Tools vignette](https://github.com/Deniz-Koseoglu/Data-Analysis-Tools/blob/master/Example%20Use%20Vignette.md) for example usage and output of this and other functions.

#### Values
Returns a list of length 3. The **first element** contains a `data.frame` of results for each pair of time series specified in `query_subset` and `ref_subset`. These include row numbers of matched `ref_data` and `query_data` (first 2 columns), x-axis values (columns 3 and 4), y-axis values (columns 5 and 6), alignment pair-specific and cumulative costs (columns 7 and 8), and the number of steps taken for each alignment pair (final column). The **second element** contains `dtw` objects for each `query_subset` - `ref_subset` pair. Finally, the **third element** is a list of `ggplot2` plots produced (**NOTE** that plots using R base graphics are *not* included).
