# 1) Principal Component Analysis (PCA, with supplementary variables and/or observations supported!)
Extended_PCA <- function(data, vars, obs="default",
                         plot_grphcs=list(c("arrow", "text"), c("point")), plot_colvar="red",
                         plot_cols, rank_vars=TRUE, add_opts=c(TRUE, FALSE, TRUE),
                         ind_groups=c("HC", "complete"), clust_k="auto", clust_dist="euclidean",
                         clust_which="none", ellipses=TRUE, pc_keep = 5,
                         export_path=getwd(), height=10, width=10, point_size=12 
) {
  
  #Install and load the required libraries 
  list.of.packages <- c("data.table", "xlsx", "FactoMineR", "factoextra", "ggplot2", "RColorBrewer", "NbClust", "dbscan", "fpc")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  invisible(lapply(c(list.of.packages, "parallel"), require, character.only = TRUE))
  
  library(data.table)
  library(xlsx)
  library(FactoMineR)
  library(factoextra)
  library(ggplot2)
  library(RColorBrewer)
  library(NbClust)
  library(dbscan)
  library(fpc)
  
  #Running checks...
  if(!all(clust_which %in% c("var","ind", "none"))) {
    print("Incorrect values provided for clust_which... Defaulting to variables only!")
    clust_which <- "var"
  }
  
  #Loading the data (no special treatment of NULL values since they are automatically deleted)
  if(is.character(data)) {
    orig_data <- main_data <- fread(data, data.table=FALSE)
    print(paste0("Data successfully loaded from: ", data, "!"))
  } else if (is.object(data)) {
    orig_data <- main_data <- data
    print("Data successfully loaded from the R environment!")
  }
  
  #Preparing variables (primary and supplementary)
  if(is.list(vars) & length(vars)==2) {
    
    pr_vars <- vars[[1]]
    sec_vars <- vars[[2]]
    
    if(any(vars[[1]] %in% vars[[2]])) {
      print(paste("Removing", length(intersect(vars[[1]], vars[[2]])), "repeated elements from the supplementary variables!"))
      sec_vars <- sec_vars[!sec_vars %in% intersect(vars[[1]], vars[[2]])]
      print("Done!")
    }
    
    all_vars <- c(pr_vars, sec_vars)
    
  } else if(is.character(vars)) {
    print("No secondary variables provided! Using only primary variables...")
    pr_vars <- vars
    sec_vars <- NULL
    all_vars <- c(pr_vars, sec_vars)
  }
  
  #Preparing observations (primary and supplementary)
  if(obs=="default") {
    
    pr_obs <- 1:nrow(main_data)
    sec_obs <- NULL
    all_obs <- c(pr_obs, sec_obs)
    
  } else if(is.list(obs) & length(obs)==2) {
    
    pr_obs <- obs[[1]]
    sec_obs <- obs[[2]]
    
    if(any(obs[[1]] %in% obs[[2]])) {
      print(paste("Removing", length(intersect(obs[[1]], obs[[2]])), "repeated elements from the supplementary observations/individuals!"))
      sec_obs <- sec_obs[!sec_obs %in% intersect(obs[[1]], obs[[2]])]
      print("Done!")
    }
    all_obs <- c(pr_obs, sec_obs)
    
  } else if(is.numeric(obs) & !is.list(obs)) {
    print("No secondary observations/individuals provided! Using only primary observations...")
    pr_obs <- obs
    sec_obs <- NULL
    all_obs <- c(pr_obs, sec_obs)
  }
  
  #Filtering main_data
  main_data <- main_data[all_obs, all_vars]
  
  if(any(is.na(main_data))) {
    print("Removing observations with missing data...")
    main_data <- na.omit(main_data)
  }
  
  #Ranking the variables to imitate Spearman correlation-based PCA!
  if(rank_vars==TRUE) {
    print("Deriving ranks for variables since rank_vars=TRUE...")
    rank_list <- as.data.frame(sapply(colnames(main_data), function(x) rank(main_data[,x], ties.method = "average")))
    main_data <- do.call(cbind, rank_list)
    print("Done!")
  }
  
  #Carrying out PCA and aggregating/exporting results (as an .xlsx file)
  if(!is.null(sec_vars)) {
    pca_res <- PCA(main_data, scale.unit = add_opts[1],
                   ncp = pc_keep, ind.sup = sec_obs,
                   quanti.sup = which(colnames(main_data) %in% sec_vars),
                   quali.sup = NULL,
                   graph = FALSE)
  } else if(is.null(sec_vars)) {
    pca_res <- PCA(main_data, scale.unit = add_opts[1],
                   ncp = pc_keep, ind.sup = sec_obs,
                   quanti.sup = NULL,
                   quali.sup = NULL,
                   graph = FALSE)
  }
  
  PCA_Res <- createWorkbook()
  
  if(nrow(pca_res[["eig"]])>pca_res[["call"]][["ncp"]]) {
    eig_rows <- pca_res[["call"]][["ncp"]]
  } else if(nrow(pca_res[["eig"]])<=pca_res[["call"]][["ncp"]]) {
    eig_rows <- nrow(pca_res[["eig"]])
  }
  
  eig_vals <- as.data.frame(pca_res[["eig"]][1:eig_rows,])
  rownames(eig_vals) <- pca_rownames <- paste0("PC", 1:eig_rows)
  
  resnames <- c("coord.pr", "coord.sec", "cos2.pr", "cos2.sec", "contrib.pr", "contrib.sec")
  x_list <- c("coord", "contrib", "cos2")
  x_list_s <- x_list[c(1,3)]
  
  reslist_prvars <- lapply(x_list, function(x) as.data.frame(t(as.data.frame(pca_res[["var"]][x]))))
  for(i in seq_along(reslist_prvars)) {
    colnames(reslist_prvars[[i]]) <- paste0(colnames(reslist_prvars[[i]]), ".", x_list[i])
    rownames(reslist_prvars[[i]]) <- pca_rownames
  }
  reslist_prvars <- do.call(cbind, reslist_prvars)
  sheetname <- "Primary Variables"
  sheet <- createSheet(PCA_Res, paste(sheetname))
  addDataFrame(reslist_prvars, sheet=sheet, startColumn=1, row.names=TRUE)
  
  
  if(any(names(pca_res) %in% "quanti.sup")) {
    reslist_secvars <- lapply(x_list[c(1,3)], function(x) as.data.frame(t(as.data.frame(pca_res[["quanti.sup"]][x]))))
    
    for(i in seq_along(reslist_secvars)) {
      colnames(reslist_secvars[[i]]) <- paste0(colnames(reslist_secvars[[i]]), ".", x_list_s[i])
      rownames(reslist_secvars[[i]]) <- pca_rownames
    }
    reslist_secvars <- do.call(cbind, reslist_secvars)
    sheetname <- "Secondary Variables"
    sheet <- createSheet(PCA_Res, paste(sheetname))
    addDataFrame(reslist_secvars, sheet=sheet, startColumn=1, row.names=TRUE)
  }
  
  reslist_prinds <- lapply(x_list, function(x) as.data.frame(pca_res[["ind"]][x]))
  #names(reslist_prinds) <- x_list
  for(i in seq_along(reslist_prinds)) {
    colnames(reslist_prinds[[i]]) <- paste0(pca_rownames, ".", x_list[i])
  }
  reslist_prinds <- do.call(cbind, reslist_prinds)
  sheetname <- "Primary Observations"
  sheet <- createSheet(PCA_Res, paste(sheetname))
  addDataFrame(reslist_prinds, sheet=sheet, startColumn=1, row.names=TRUE)
  
  if(any(names(pca_res) %in% "ind.sup")) {
    reslist_secinds <- lapply(x_list[c(1,3)], function(x) as.data.frame(pca_res[["ind.sup"]][x]))
    #names(reslist_secinds) <- x_list[c(1,3)]
    for(i in seq_along(reslist_secinds)) {
      colnames(reslist_secinds[[i]]) <- paste0(pca_rownames, ".", x_list[i])
    }
    reslist_secinds <- do.call(cbind, reslist_secinds)
    sheetname <- "Secondary Observations"
    sheet <- createSheet(PCA_Res, paste(sheetname))
    addDataFrame(reslist_secinds, sheet=sheet, startColumn=1, row.names=TRUE)
  }
  
  #Exporting all results as an .xlsx file
  saveWorkbook(PCA_Res, paste0(export_path, "/EasyPCA_Results_", format(Sys.time(), "%Y-%m-%d %Hhr %Mmin %Ssec"), ".xlsx"))
  
  #Creating plots
  pca_plotlist <- list()
  
  #1) Scree plot of eigen values
  pca_plotlist[["eigen_scree"]] <- fviz_eig(pca_res, addlabels = TRUE, ylim = c(0, 100))
  
  if(length(ind_groups)==1 & any(ind_groups[[1]] %in% c("HC", "kmeans"))) {
    
    if(ind_groups[[1]]=="HC") {
      print("No linkage method for HC was provided! Defaulting to complete linkage...")
      ind_groups[[2]] <- c("complete")
      
    } else if(ind_groups[[1]]=="kmeans") {
      
      ind_groups[[2]] <- "kmeans"
      
    } 
    
  } else if(length(ind_groups)<3 & ind_groups[[1]]=="dbscan") {
    
    if(length(ind_groups)==2 & !any(ind_groups[[2]]=="auto") & !is.numeric(ind_groups[[2]])) {
      print("Incorrect EPS value provided for dbscan! Switching to default...")
      ind_groups[[2]] <- "auto"
      ind_groups[[3]] <- 5
    }
    
  }
  
  #Running Final Checks
  if(ind_groups[[1]]=="HC") {
    
    if(!any(ind_groups[[2]] %in% c("single", "average", "complete", "ward.D", "ward.D2", "mcquitty"))) {
      print("Incorrect HC linkage method provided! Defaulting to complete linkage...")
      ind_groups[[2]] <- "complete"
    }
    
  } else if(ind_groups[[1]]=="kmeans") {
    
    if(ind_groups[[2]]!="kmeans") {
      ind_groups[[2]] <- "kmeans"
    }
    
  } else if(ind_groups[[1]]=="dbscan") {
    
    if(!any(ind_groups[[2]]=="auto") & !is.numeric(as.numeric(ind_groups[[2]]))) {
      print("Incorrect EPS value provided for DBSCAN! Swithing to automatic determination...")
      ind_groups[[2]] <- "auto"
    }
    if(!is.numeric(as.numeric(ind_groups[[3]]))) {
      print("Only integer values are supported for DBSCAN MinPts! Defaulting to 5...")
      ind_groups[[3]] <- 5
    }
    
  }
  
  #Assessing clustering tendency of data (Hopkins statistic and HeatMap)
  hop_stat <- c()
  clust_whichlab <- c()
  if(any(clust_which %in% "ind")) { clust_whichlab[which(clust_which=="ind")] <- "Individuals" }
  if(any(clust_which %in% "var")) { clust_whichlab[which(clust_which=="var")] <- "Variables" }
  
  for(i in seq_along(clust_which)) {
    hop_stat[clust_whichlab[i]] <- factoextra::get_clust_tendency(as.data.frame(pca_res[[clust_which[i]]]["coord"]), 
                                                                  n = nrow(as.data.frame(pca_res[[clust_which[i]]]["coord"]))-1, 
                                                                  graph = FALSE)[["hopkins_stat"]]
    
    pca_plotlist[[paste0("ClusTend_HeatMap_", clust_whichlab[i])]] <- fviz_dist(dist(as.data.frame(pca_res[[clust_which[i]]]["coord"]), method=clust_dist), 
                                                                                show_labels = FALSE) + labs(title = paste(clust_whichlab[i], "clustering tendency"), caption = paste0("Hopkins Statistic: ", hop_stat[clust_whichlab[i]]))
  }
  
  #Optimal number of clusters
  if(length(clust_which)>1 & length(clust_k)==1) {
    clust_k <- rep(clust_k, length(clust_which))
  }
  
  if(length(clust_which)>1) {
    if(ind_groups[[1]]=="dbscan") {
      if(length(ind_groups[[2]])==1) {
        ind_groups[[2]] <- rep(ind_groups[[2]], length(clust_which))
      }
    }
    if(ind_groups[[1]]=="dbscan") {
      if(length(ind_groups[[3]])==1) {
        ind_groups[[3]] <- rep(ind_groups[[2]], length(clust_which))
      }
    }
  }
  
  opt_cnum <- c()
  
  if(!any(clust_which %in% "none")) {
    
    cluster_res <- list()
    
    if(ind_groups[[1]]!=FALSE & !any(ind_groups[[1]] %in% c("HC", "kmeans", "dbscan", "fcm"))) {
      print("Custom column selected for colouring observations. Processing...")
      opt_cnum["ind"] <- as.factor(orig_data[,ind_groups[[1]]])
      
    } else if(ind_groups[[1]]!=FALSE & all(ind_groups[[1]] %in% c("HC", "kmeans"))) {
      
      #optname <- paste0("Optimal Cluster Number_", ind_groups[2])
      optimal_clusters <- list()
      
      for(i in seq_along(clust_which)) {
        if(clust_k[i]=="auto") {
          print(paste("Deriving optimal cluster number for method:", ind_groups[[1]], "using NbClust..."))
          
          print(paste0("Working with ", clust_which[i], "..."))
          pdf(file=NULL)
          optimal_clusters[[clust_which[i]]] <- try(NbClust(as.data.frame(pca_res[[clust_which[i]]]["coord"]), max.nc = pmax(2, pmin(nrow(as.data.frame(pca_res[[clust_which[i]]]["coord"])), 16))-1, 
                                                            distance = clust_dist, method = ind_groups[[2]], index="all"))
          dev.off()
          if(inherits(optimal_clusters[[clust_which[i]]], "try-error")) {
            
            print("At least one NbClust optimal cluster test failed! Repeating procedure using individual indices...")
            pdf(file=NULL)
            
            index_list <- c("kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", "friedman", "rubin", "cindex", "db", "silhouette", "duda",
                            "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", "gamma", "gplus", "tau", "dunn", "hubert",
                            "sdindex", "dindex", "sdbw")
            
            optimal_clusters[[clust_which[i]]] <- sapply(index_list, function(x) as.numeric(try(NbClust(as.data.frame(pca_res[[clust_which[i]]]["coord"]), max.nc = pmax(2, pmin(nrow(as.data.frame(pca_res[[clust_which[i]]]["coord"])), 16))-1, 
                                                                                                        distance = clust_dist, method = ind_groups[[2]], index= x)[["Best.nc"]]["Number_clusters"])))
            dev.off()
            
            failed_tests <- unique(c(which(sapply(names(optimal_clusters[[clust_which[i]]]), function(x) length(optimal_clusters[[clust_which[i]]][[x]]))==0), which(is.na(optimal_clusters[[clust_which[i]]]))))
            
            if(length(failed_tests)==length(index_list)) {
              
              print("All NbClust indices failed to converge! Defaulting to 3 clusters...")
              opt_cnum[clust_which[i]] <- 3
              
            } else if(length(failed_tests)!=length(index_list)) {
              
              print(paste0("A total of ", length(failed_tests), " NbClust tests failed to converge! Using the remaining ", length(index_list)-length(failed_tests), "..."))
              
              freq_table <- table(unlist(optimal_clusters[[clust_which[i]]])[-failed_tests])
              opt_cnum[clust_which[i]] <- as.numeric(names(freq_table)[freq_table==max(freq_table)])[1]
              
            }
            
          } else if (!inherits(optimal_clusters[[clust_which[i]]], "try-error")) {
            opt_cnum[clust_which[i]] <- as.numeric(names(which.max(summary(as.factor(optimal_clusters[[clust_which[i]]][["Best.nc"]][1,])))))
            print(paste("Determined an optimal number of clusters to be", opt_cnum, "by majority rule!"))
          }
        } else if(is.numeric(clust_k[i])) {
          print("Using manually provided number of clusters...")
          #opt_cnum <- rep(clust_k, length(clust_which))
          opt_cnum[i] <- clust_k[i]
          names(opt_cnum)[i] <- clust_which[i]
        }
        
      } 
      
      if(ind_groups[[1]]=="HC") {
        for(i in seq_along(clust_which)) {
          dist_obj <- dist(as.data.frame(pca_res[[clust_which[i]]]["coord"]), method=clust_dist)
          cluster_obj <- hclust(dist_obj, method=ind_groups[[2]])
          cluster_res[[clust_which[i]]] <- as.factor(cutree(cluster_obj, k=opt_cnum[i]))
        }
      } else if(ind_groups[[1]]=="kmeans") {
        for(i in seq_along(clust_which)) {
          cluster_res[[clust_which[i]]] <- as.factor(kmeans(pca_res[[clust_which[i]]][["coord"]],
                                                            centers=opt_cnum[i],
                                                            nstart=25,
                                                            iter.max=500)[["cluster"]])
        }
      }
      
    } else if(ind_groups[[1]]!=FALSE & ind_groups[[1]]=="dbscan") {
      
      print("NOTE: Cluster number cannot be pre-defined for DBSCAN! Argument clust_k will be ignored...")
      eps_val <- c()
      
      for(i in seq_along(clust_which)) {
        if(ind_groups[[2]][i]=="auto") {
          print(paste0("Attempting to determine appropriate EPS automatically for ", clust_whichlab[i], "..."))
          db_kNN <- dbscan::kNNdist(as.data.frame(pca_res[[clust_which[i]]]["coord"]), k = ind_groups[[3]][i], all = FALSE)
          db_kNN <- db_kNN[order(db_kNN)]
          
          d1 <- diff(db_kNN)/diff(seq_along(db_kNN))
          d2 <- diff(d1)/diff(seq_along(db_kNN)[-1])
          
          eps_val[i] <- max(db_kNN[which.max(d2)])
          
        } else if(ind_groups[[2]][i]!="auto") {
          print(paste0("Using provided EPS value for ", clust_whichlab[i], "..."))  
          eps_val[i] <- ind_groups[[2]][i]
        }
        
        cluster_res[[clust_which[i]]] <- as.factor(fpc::dbscan(as.data.frame(pca_res[[clust_which[i]]]["coord"]),
                                                               eps = eps_val[i],
                                                               MinPts = ind_groups[[3]][i],
                                                               method = "raw")[["cluster"]])
        
      }
      
    } else if(ind_groups[[1]]!=FALSE & ind_groups[[1]]=="fcm") {
      
      if(any(clust_k=="auto")) {
        print("Automatic cluster detection is NOT supported for FCM! Defaulting to 3 clusters...")
        clust_k[which(clust_k=="auto")] <- 3
        clust_k <- as.numeric(clust_k)
      }
      
      print("Using manually provided number of clusters...")
      
      clust_probs <- list()
      
      for(i in seq_along(clust_which)) {
        dist_obj <- dist(as.data.frame(pca_res[[clust_which[i]]]["coord"]), method=clust_dist)
        clust_obj <- fanny(dist_obj, k = clust_k[i], diss = TRUE)
        cluster_res[[clust_which[i]]] <- as.factor(clust_obj[["clustering"]])
        clust_probs[[clust_which[i]]] <- as.data.frame(clust_obj[["membership"]])
        colnames(clust_probs[[clust_which[i]]]) <- paste0("Cluster_", 1:ncol(clust_probs[[clust_which[i]]]))
      }
    }
  } else cluster_res <- list()
  
  #Setting up some PCA args
  #Setting up geoms
  geom_var <- plot_grphcs[[1]]
  geom_ind <- plot_grphcs[[2]]
  
  #Setting up other plot options
  col_ind <- "black"
  col_ind_sup <- plot_cols[[4]]
  col_var_sup <- plot_cols[[3]]
  mean_point <- add_opts[2]
  repel_text <- add_opts[3]
  #gradient_cols <- plot_cols[[1]]
  #palette_cols <- plot_cols[[2]]
  point_indshape <- 21
  point_varshape <- 1
  #label = "all"
  
  #Setting up main argument-dependent PCA plot options
  if(any(names(cluster_res) %in% "var")) {
    
    col_var <- cluster_res[["var"]]
    
    pca_plotlist[["PCA_Variables"]] <- fviz_pca_var(pca_res, geom.var = geom_var,
                                                    repel = repel_text,
                                                    col.var = col_var,
                                                    col.quanti.sup = col_var_sup,
                                                    palette = plot_cols[[1]],
                                                    col.circle = "black")
    
  } else if(!any(names(cluster_res) %in% "var")) {
    
    col_var <- plot_colvar
    
    pca_plotlist[["PCA_Variables"]] <- fviz_pca_var(pca_res, geom.var = geom_var,
                                                    repel = repel_text,
                                                    col.var = col_var,
                                                    col.quanti.sup = col_var_sup,
                                                    gradient.cols = plot_cols[[1]],
                                                    col.circle = "black")
    
  }
  
  if(any(names(cluster_res) %in% "ind")) {
    
    #habil <- cluster_res[["ind"]]
    col_indfill <- cluster_res[["ind"]]
    add_ell_switch <- ellipses
    
    pca_plotlist[["PCA_Individuals"]] <- fviz_pca_ind(pca_res, geom.ind = geom_ind,
                                                      repel = repel_text,
                                                      #habillage = habil,
                                                      palette = plot_cols[[2]],
                                                      mean.point = mean_point,
                                                      addEllipses = add_ell_switch,
                                                      #ellipse.type = add_ellipses,
                                                      pointshape = point_indshape,
                                                      col.ind = col_ind,
                                                      col.ind.sup = col_ind,
                                                      fill.ind = col_indfill,
                                                      fill.ind.sup = col_ind_sup)
    
  } else if(!any(names(cluster_res) %in% "ind")) {
    
    #habil <- "none"
    col_indfill <- plot_cols[[5]]
    add_ell_switch <- FALSE
    
    pca_plotlist[["PCA_Individuals"]] <- fviz_pca_ind(pca_res, geom.ind = geom_ind,
                                                      repel = repel_text,
                                                      #habillage = habil,
                                                      #gradient.cols = plot_cols[[2]],
                                                      mean.point = mean_point,
                                                      addEllipses = add_ell_switch,
                                                      #ellipse.type = add_ellipses,
                                                      pointshape = point_indshape,
                                                      col.ind = col_ind,
                                                      col.ind.sup = col_ind,
                                                      fill.ind = col_indfill,
                                                      fill.ind.sup = col_ind_sup)
    
  }
  
  if(any(clust_which %in% "var")) {
    pca_plotlist[["PCA_Biplot"]] <- fviz_pca_biplot(pca_res, geom.var = geom_var,
                                                    geom.ind = geom_ind,
                                                    repel = repel_text,
                                                    mean.point.ind = mean_point,
                                                    addEllipses = FALSE,
                                                    pointshape = point_indshape,
                                                    col.ind = col_ind,
                                                    col.ind.sup = col_ind,
                                                    fill.ind = col_indfill,
                                                    fill.ind.sup = col_ind_sup,
                                                    col.var = col_var,
                                                    col.quanti.sup = col_var_sup) +
      ggpubr::fill_palette(plot_cols[[2]]) +      # Indiviual fill color
      ggpubr::color_palette(plot_cols[[1]])      # Variable colors
  }
  
  pca_plotlist[["Individuals_Contribution_PCs1_and_2"]] <- fviz_contrib(pca_res, choice = "ind", axes = 1:2)
  pca_plotlist[["Variables_Contribution_PC1"]] <- fviz_contrib(pca_res, choice = "var", axes = 1)
  pca_plotlist[["Variables_Contribution_PC2"]] <- fviz_contrib(pca_res, choice = "var", axes = 2)
  pca_plotlist[["Variables_Contribution_PCs1_and_2"]] <- fviz_contrib(pca_res, choice = "var", axes = 1:2)
  
  print("Exporting plots as a .PDF...")
  pdf(paste0(export_path, "/EasyPCA_Plots_", format(Sys.time(), "%Y-%m-%d %Hhr %Mmin %Ssec"), ".pdf"),
      width = width,
      height = height,
      pointsize = point_size)
  print(pca_plotlist)
  dev.off()
  
  if(any(clust_which %in% c("var", "ind"))) {
    print("Exporting plot data as and .XLSX workbook...")
    PCA_Plotres <- createWorkbook()
    
    if(any(clust_which %in% "var")) {
      
      sheet <- createSheet(PCA_Plotres, "Primary Variables")
      
      if(ind_groups[[1]]=="fcm") {
        clust_totalres_var <- cbind.data.frame(pca_plotlist[["PCA_Variables"]][["data"]], clust_probs[["var"]])
      } else if(ind_groups[[1]]!="fcm") {
        clust_totalres_var <- as.data.frame(pca_plotlist[["PCA_Variables"]][["data"]]) 
      }
      
      addDataFrame(clust_totalres_var, sheet=sheet, startColumn=1, row.names=TRUE)
    }
    
    if(any(clust_which %in% "ind")) {
      
      sheet <- createSheet(PCA_Plotres, "Primary Individuals")
      
      if(ind_groups[[1]]=="fcm") {
        clust_totalres_ind <- cbind.data.frame(pca_plotlist[["PCA_Individuals"]][["data"]], clust_probs[["ind"]])
      } else if(ind_groups[[1]]!="fcm") {
        clust_totalres_ind <- as.data.frame(pca_plotlist[["PCA_Individuals"]][["data"]]) 
      }
      
      addDataFrame(clust_totalres_ind, sheet=sheet, startColumn=1, row.names=TRUE)
    }
    
    #Exporting all results as an .xlsx file
    saveWorkbook(PCA_Plotres, paste0(export_path, "/EasyPCA_Plot_Cluster_Data_", format(Sys.time(), "%Y-%m-%d %Hhr %Mmin %Ssec"), ".xlsx"))
  }
  
  print("Processing complete!")
  #print(pca_plotlist)
  
  if(ind_groups[[1]]=="fcm") {
    pca_fres <- list(pca_res, clust_probs)
  } else pca_fres <- pca_res
  
  return(list(Results=pca_fres, Plots=pca_plotlist))
}

# 2) Hierarchical Clustering (by variable or observation!)
Extended_HC <- function(data, cluster_by="observations", var_subset, var_labels=var_subset, dist_measure="euclidean",
                        agglomeration_methods=c("single", "average", "complete", "ward.D", "ward.D2", "mcquitty"),
                        agglomeration_labels=c("Single", "Average", "Complete", "ward.D", "ward.D2", "McQuitty"),
                        preproc_method="zscore", cvi_range=c(2,10), cluster_legend=TRUE, cluster_labels=NULL, k=3, cor_method="cophenetic", 
                        horizont=TRUE, dataset_label="Sample", cluster_colours=c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2")),
                        gg_lab_size=0.9, base_lab_size=0.8, gg_lwd=0.8, base_lwd=1.2, draw_rect=TRUE, coef_method="spearman",
                        FM_test=NULL, FM_lim=c(0.7,1), export_plots="pdf", width=5, height=5, point_size=10, dpi=500,
                        export_results=TRUE, export_path=getwd(), dend_mar=c(4,4)) {
  
  #Install and load the required libraries 
  list.of.packages <- c("biotools", "asbio", "readr", "lattice", "corrplot", "fmsb", "latticeExtra", "RColorBrewer",
                        "reshape2", "gridExtra", "car", "psych", "data.table", "gridGraphics", "ggcorrplot", "corrgram",
                        "Hmisc", "xlsx", "factoextra", "NbClust", "dendextend", "fpc", "dtw", "dtwclust", "TSclust", 
                        "rlist", "pryr", "GGally", "ggpmisc", "lmtest", "cpm", "ecp", "changepoint", "magrittr", "scales")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  
  library(biotools)
  library(asbio)
  library(readr)
  library(lattice)
  library(corrplot)
  library(fmsb)
  library(latticeExtra)
  library(RColorBrewer)
  library(reshape2)
  library(gridExtra)
  library(car)
  library(psych)
  library(data.table)
  library(gridGraphics)
  library(ggcorrplot)
  library(corrgram)
  library(Hmisc)
  library(xlsx)
  library(factoextra)
  library(NbClust)
  library(dendextend)
  library(fpc)
  library(dtw)
  library(dtwclust)
  library(TSclust)
  library(rlist)
  library(pryr)
  library(GGally)
  library(ggpmisc)
  library(lmtest)
  library(cpm)
  library(ecp)
  library(changepoint)
  library(magrittr)
  library(scales)
  
  if(cluster_by=="observations" & dist_measure %in% c("DTW", "GAK")) {
    stop("DTW and GAK distances are not allowed when clustering observations! Please use cluster_by = \"variables\" instead.")
  }
  
  if(!is.null(FM_test) & cluster_by!="observations") {
    warning("Fawlkes-Mallows testing cannot be carried out when clustering by variables! Skipping.")
  }
  
  if(cluster_by=="observations") {
    Data_List <- as.data.frame(data[,var_subset])
    Data_List <- na.omit(Data_List)
    Data_List <- scale(Data_List)
  } else if(cluster_by=="variables") {
    Data_List <- as.data.frame(data[,var_subset])
    Data_List <- na.omit(Data_List)
    Data_List <- as.data.frame(scale(Data_List))
    Data_List <- as.list(Data_List)
  }
  
  HC_clusters <- list()
  
  if(cluster_by=="observations") {
    HC_clusters[["Clustering Tendency Heatmap"]] <- fviz_dist(dist(as.data.frame(Data_List), method=dist_measure), show_labels = FALSE)+
      labs(title = paste(dataset_label, " clustering tendency", sep=""))
    
    GS14_Hopkins <- get_clust_tendency(as.data.frame(Data_List), n = nrow(as.data.frame(Data_List))-1, graph = FALSE)
    GS14_Hopkins$hopkins_stat
    
    HC_clusters[["Elbow Method"]] <- fviz_nbclust(as.data.frame(Data_List), hcut, diss=proxy::dist(as.data.frame(Data_List), method = dist_measure), method = "wss") +
      geom_vline(xintercept=3, linetype = 2) +
      labs(title=paste(dataset_label, " optimal cluster number", sep=""), subtitle = "Elbow method")
    
    HC_clusters[["Average Silhouette"]]<- fviz_nbclust(as.data.frame(Data_List), hcut, diss=proxy::dist(as.data.frame(Data_List), method = dist_measure), method = "silhouette")+
      labs(title=paste(dataset_label, " optimal cluster number", sep=""), subtitle = "Silhouette method")
    
    set.seed(123)
    HC_clusters[["Gap Statistic"]] <- fviz_nbclust(as.data.frame(Data_List), hcut, diss=proxy::dist(x, method = dist_measure), nstart = 25,  method = "gap_stat", nboot = 500)+
      labs(title=paste(dataset_label, " optimal cluster number", sep=""), subtitle = "Gap statistic method")
    
    GS14_optimal_clusters <- list()
    
    for(i in 1:length(agglomeration_methods)) {
      optname <- paste("Optimal Cluster Number_", agglomeration_labels[i], " Linkage", sep="")
      GS14_optimal_clusters[[optname]] <- try(NbClust(as.data.frame(Data_List), distance = dist_measure, min.nc = cvi_range[1],
                                                      max.nc = cvi_range[2], method = agglomeration_methods[i]))
      if(inherits(GS14_optimal_clusters[[optname]], "try-error")) {
        print(paste0("NbClust optimal cluster testing failed to converge for ", agglomeration_labels[i], "! The associated graph will not be plotted..."))
      } else if (!inherits(GS14_optimal_clusters[[optname]], "try-error")) {
        HC_clusters[[paste("30-rule testing_", agglomeration_labels[i], " Linkage", sep="")]] <- fviz_nbclust(GS14_optimal_clusters[[optname]]) + 
          labs (subtitle = paste("30-rule testing_", agglomeration_labels[i], " Linkage", sep="")) + 
          theme(panel.grid.major.y = element_line(size=0.5, color="grey75", linetype=3), plot.subtitle = element_text())
        
      }
    }
    
    if(export_results==TRUE & cluster_by=="observations") {
      CVI_Results <- createWorkbook()
      for(i in seq_along(GS14_optimal_clusters)) {
        sheetname <- paste(names(GS14_optimal_clusters[i]))
        sheet <- createSheet(CVI_Results, paste(sheetname))
        addDataFrame(GS14_optimal_clusters[[i]]$Best.nc, sheet=sheet, startColumn=1, row.names=FALSE)
      }
      saveWorkbook(CVI_Results, paste0(export_path, "/CVI Results (by Observations)_", format(Sys.time(), "%Y-%m-%d %Hhr %Mmin %Ssec"), ".xlsx"))
    }
    
  }
  
  Data_HC_clusters <- list()
  
  if(!is.null(preproc_method)) {
    for(i in 1:length(agglomeration_methods)) {
      HC_listname <- paste(agglomeration_labels[i], sep="")
      Data_HC_clusters[[HC_listname]] <- tsclust(Data_List, type="hierarchical",
                                                 k=k, seed=42, preproc=eval(parse(text=preproc_method)),
                                                 distance=dist_measure, control=hierarchical_control(method=agglomeration_methods[i]),
                                                 error.check=TRUE) # args=tsclust_args(dist=list(sigma=sigma, window.size=18L))
    }
  } else if(is.null(preproc_method)) {
    for(i in 1:length(agglomeration_methods)) {
      HC_listname <- paste(agglomeration_labels[i], sep="")
      Data_HC_clusters[[HC_listname]] <- tsclust(Data_List, type="hierarchical",
                                                 k=k, seed=42,
                                                 distance=dist_measure, control=hierarchical_control(method=agglomeration_methods[i]),
                                                 error.check=TRUE) # args=tsclust_args(dist=list(sigma=sigma, window.size=18L))
    }
  }
  
  #Deriving cluster memberships from AHC trees
  cluster_memberships <- lapply(Data_HC_clusters, function(x) as.data.frame(stats::cutree(x, k = k)))
  cluster_memberships <- do.call(cbind.data.frame, cluster_memberships)
  colnames(cluster_memberships) <- names(Data_HC_clusters)
  
  if(export_results==TRUE) {
    fwrite(cluster_memberships, file = paste0(export_path, "/Cluster_MShips_", format(Sys.time(), "%Y-%m-%d %Hhr %Mmin %Ssec"), ".csv"), na = NA)
  }
  
  #if(cluster_by=="variables") {
  #  cvi_data <- list()
  
  #  for(i in 1:length(agglomeration_methods)) {
  #    HC_listname <- paste(agglomeration_labels[i], sep="")
  #    cvi_result <- tsclust(Data_List, type="hierarchical",
  #                                               k=cvi_low:cvi_high, seed=42, #preproc=zscore
  #                                               distance=dist_measure, control=hierarchical_control(method=agglomeration_methods[i]),
  #                                               error.check=TRUE) # args=tsclust_args(dist=list(sigma=sigma, window.size=18L))
  #    names(cvi_result) <- paste0("k_",cvi_low:cvi_high)
  #    cvi_data[[HC_listname]] <- sapply(cvi_result, cvi, type="internal")
  #    }
  #}
  
  #if(export_results==TRUE & cluster_by=="variables") {
  #  write.csv(cvi_results, "Cluster Evaluation Indices (by variables).csv", row.names=FALSE)
  #}
  
  for(i in 1:length(Data_HC_clusters)) {
    
    Data_HC_clusters[[i]][["labels"]] <- var_labels
    
    GG_HC_listname <- paste("GGPlot_",agglomeration_labels[i], sep="")
    HC_clusters[[GG_HC_listname]] <- fviz_dend(Data_HC_clusters[[i]], k = k, # Cut in four groups
                                               cex = gg_lab_size, # label size
                                               k_colors = cluster_colours,
                                               color_labels_by_k = TRUE, # color labels by groups
                                               rect = draw_rect, # Add rectangle around groups
                                               type="rectangle",
                                               horiz=horizont,
                                               ylab=paste("Dissimilarity (", dist_measure,")", sep=""),
                                               xlab="",
                                               main=paste(dataset_label, " clusters (", agglomeration_labels[i]," linkage)", sep=""),
                                               lwd=gg_lwd,
                                               ggtheme=theme(axis.line.x=element_line(colour="black", size=0.8),
                                                             axis.ticks.x=element_line(colour="black", size=0.8),
                                                             axis.text = element_text(size=11),
                                                             axis.title = element_text(size=11),
                                                             panel.background=element_blank(),
                                                             panel.grid=element_blank(),
                                                             legend.justification= c(1, 0), 
                                                             legend.position = c(1, 0.5),
                                                             panel.border=element_rect(colour="black", size=0.5, fill=NA)))
  }
  
  # Another method to plot a dendrogram using the dendextend package
  HC_dendextend <- list()
  
  for(i in 1:length(Data_HC_clusters)) {
    dendextend_name <- paste("Base_",agglomeration_labels[i], sep="")
    HC_dendextend[[dendextend_name]] <- as.dendrogram(Data_HC_clusters[[i]])
    
    if(horizont==FALSE) {
      #pryr_plot1 %<a-% {
      par(mgp=c(1.6,0.4,0),mar=c(dend_mar[1],4,2,dend_mar[2]), mfrow=c(1,1))
      HC_dendextend[[dendextend_name]] %>%
        color_branches(.,k=k, col=cluster_colours) %>%
        set("labels_col", value=cluster_colours, k=k) %>%
        set("labels_cex", base_lab_size) %>%
        set("branches_lwd", base_lwd) %>% plot(horiz=horizont, yaxt="n", ylab=paste("Dissimilarity (", dist_measure,")", sep=""),
                                               main=paste(dataset_label, " clusters (", agglomeration_labels[i], " linkage)", sep=""))
      
      axis(2, cex.axis=0.9, tck=-0.013) 
      
      if(draw_rect==TRUE) {
        HC_dendextend[[dendextend_name]] %>% rect.dendrogram(k=k, horiz = horizont,
                                                             border = "grey70", lty = 5, lwd = base_lwd)
      }
      
      HC_dendextend[[dendextend_name]] %<>% set("labels_col", value = cluster_colours, k= k)
      
      if(cluster_legend==TRUE) {
        # For the cluster legend...
        categories <- cutree(HC_dendextend[[dendextend_name]], k = k)
        colors_HC <- labels_colors(HC_dendextend[[dendextend_name]])[categories %>% sort %>% names]
        clust_labs <- colors_HC %>% unique
        clust_labs <- clust_labs[order(match(clust_labs, cluster_colours[1:length(clust_labs)]))]
        clust_leg <- unique(categories)
        
        if(!is.null(cluster_labels)) {
          clust_leg <- cluster_labels
        }
        legend("topright", legend = clust_leg, fill = clust_labs, cex = 0.9, bty="n")
      }
      
      if(export_plots=="png") {
        png(filename=paste(export_path,"/Base_HC_",agglomeration_labels[i], "_", format(Sys.time(), "%Y-%m-%d %Hhr %Mmin %Ssec"), ".png", sep=""),
            type="cairo", units="in",
            width=width, height=height,
            pointsize=point_size,
            res=dpi)
        par(mgp=c(1.6,0.4,0),mar=c(dend_mar[1],4,2,dend_mar[2]), mfrow=c(1,1))
        HC_dendextend[[dendextend_name]] %>%
          color_branches(.,k=k, col=cluster_colours) %>%
          set("labels_col", value=cluster_colours, k=k) %>%
          set("labels_cex", base_lab_size) %>%
          set("branches_lwd", base_lwd) %>% plot(horiz=horizont, yaxt="n", ylab=paste("Dissimilarity (", dist_measure,")", sep=""),
                                                 main=paste(dataset_label, " clusters (", agglomeration_labels[i], " linkage)", sep=""))
        
        axis(2, cex.axis=0.9, tck=-0.013) 
        
        if(draw_rect==TRUE) {
          HC_dendextend[[dendextend_name]] %>% rect.dendrogram(k=k, horiz = horizont,
                                                               border = "grey70", lty = 5, lwd = base_lwd)
        }
        
        HC_dendextend[[dendextend_name]] %<>% set("labels_col", value = cluster_colours, k= k)
        
        if(cluster_legend==TRUE) {
          # For the cluster legend...
          categories <- cutree(HC_dendextend[[dendextend_name]], k = k)
          colors_HC <- labels_colors(HC_dendextend[[dendextend_name]])[categories %>% sort %>% names]
          clust_labs <- colors_HC %>% unique
          clust_labs <- clust_labs[order(match(clust_labs, cluster_colours[1:length(clust_labs)]))]
          clust_leg <- unique(categories)
          
          if(!is.null(cluster_labels)) {
            clust_leg <- cluster_labels
          }
          legend("topright", legend = clust_leg, fill = clust_labs, cex = 0.9, bty="n")
        }
        dev.off()
      }
      #}
      #HC_clusters[[dendextend_name]] <- pryr_plot1
    } else if(horizont==TRUE) {
      #pryr_plot1 %<a-% {
      par(mgp=c(1.6,0.4,0),mar=c(dend_mar[1],2,2,dend_mar[2]), mfrow=c(1,1))
      HC_dendextend[[dendextend_name]] %>%
        color_branches(.,k=k, col=cluster_colours) %>%
        set("labels_col", value=cluster_colours, k=k) %>%
        set("labels_cex", base_lab_size) %>%
        set("branches_lwd", base_lwd) %>% plot(horiz=horizont, yaxt="n", xaxt="n", xlab=paste("Dissimilarity (", dist_measure,")", sep=""),
                                               main=paste(dataset_label, " clusters (", agglomeration_labels[i], " linkage)", sep=""))
      
      axis(1, cex.axis=0.9, tck=-0.013) 
      
      if(draw_rect==TRUE) {
        HC_dendextend[[dendextend_name]] %>% rect.dendrogram(k=k, horiz = horizont,
                                                             border = "grey70", lty = 5, lwd = base_lwd)
      }
      
      HC_dendextend[[dendextend_name]] %<>% set("labels_col", value = cluster_colours, k = k)
      
      if(cluster_legend==TRUE) {
        # For the cluster legend...
        categories <- cutree(HC_dendextend[[dendextend_name]], k = k)
        colors_HC <- labels_colors(HC_dendextend[[dendextend_name]])[categories %>% sort %>% names]
        clust_labs <- colors_HC %>% unique
        clust_labs <- clust_labs[order(match(clust_labs, cluster_colours[1:length(clust_labs)]))]
        clust_leg <- unique(categories)
        
        if(!is.null(cluster_labels)) {
          clust_leg <- cluster_labels
        }
        legend("topleft", legend = clust_leg, fill = clust_labs, cex = 0.9, bty="n")
      }
      
      if(export_plots=="png") {
        png(filename=paste(export_path, "/Base_HC_",agglomeration_labels[i], "_", format(Sys.time(), "%Y-%m-%d %Hhr %Mmin %Ssec"), ".png", sep=""),
            type="cairo", units="in",
            width=width, height=height,
            pointsize=point_size,
            res=dpi)
        
        par(mgp=c(1.6,0.4,0),mar=c(dend_mar[1],2,2,dend_mar[2]), mfrow=c(1,1))
        HC_dendextend[[dendextend_name]] %>%
          color_branches(.,k=k, col=cluster_colours) %>%
          set("labels_col", value=cluster_colours, k=k) %>%
          set("labels_cex", base_lab_size) %>%
          set("branches_lwd", base_lwd) %>% plot(horiz=horizont, yaxt="n", xaxt="n", xlab=paste("Dissimilarity (", dist_measure,")", sep=""),
                                                 main=paste(dataset_label, " clusters (", agglomeration_labels[i], " linkage)", sep=""))
        
        axis(1, cex.axis=0.9, tck=-0.013) 
        
        if(draw_rect==TRUE) {
          HC_dendextend[[dendextend_name]] %>% rect.dendrogram(k=k, horiz = horizont,
                                                               border = "grey70", lty = 5, lwd = base_lwd)
        }
        
        HC_dendextend[[dendextend_name]] %<>% set("labels_col", value = cluster_colours, k = k)
        
        if(cluster_legend==TRUE) {
          # For the cluster legend...
          categories <- cutree(HC_dendextend[[dendextend_name]], k = k)
          colors_HC <- labels_colors(HC_dendextend[[dendextend_name]])[categories %>% sort %>% names]
          clust_labs <- colors_HC %>% unique
          clust_labs <- clust_labs[order(match(clust_labs, cluster_colours[1:length(clust_labs)]))]
          clust_leg <- unique(categories)
          
          if(!is.null(cluster_labels)) {
            clust_leg <- cluster_labels
          }
          legend("topleft", legend = clust_leg, fill = clust_labs, cex = 0.9, bty="n")
        }
        
        dev.off()
      }
      #HC_clusters[[dendextend_name]] <- pryr_plot1
      #}
    }
  }
  
  HC_cor <- cor.dendlist(as.dendlist(HC_dendextend), method=cor_method, method_coef=coef_method)
  colnames(HC_cor) <- agglomeration_labels
  rownames(HC_cor) <- as.character(c(1:length(agglomeration_labels)))
  corrplot.mixed(HC_cor, upper="number", lower="pie", tl.col="black", tl.cex=0.8, number.cex=0.8, pch.cex=1.5, title=paste("Linkage Method Comparison (", dist_measure, " distance)", sep=""), mar=c(0,0,1,0))
  
  if(export_plots=="png") {
    png(filename=paste(export_path, "/Linkage Method Correlations_", format(Sys.time(), "%Y-%m-%d %Hhr %Mmin %Ssec"), ".png", sep=""),
        type="cairo", units="in",
        width=width, height=height,
        pointsize=point_size,
        res=dpi)
    corrplot.mixed(HC_cor, upper="number", lower="pie", tl.col="black", tl.cex=0.8, number.cex=0.8, pch.cex=1.5, title=paste("Linkage Method Comparison (", dist_measure, " distance)", sep=""), mar=c(0,0,1,0))
    dev.off()
  }
  
  #pryr_plot2 %<a-% {}
  #HC_clusters[["Correlations (by Linkage Method)"]] <- pryr_plot2
  
  if(cluster_by=="observations" & !is.null(FM_test)) {
    get_ordered_clusters <- function(dend, kc) {
      cutree(dend, k = kc)[order.dendrogram(dend)]
    }
    
    HC_dendlist <- as.dendlist(HC_dendextend)
    names(HC_dendlist) <- agglomeration_labels
    dend_clusters <- lapply(HC_dendlist, FUN=get_ordered_clusters, kc = k)
    
    Data_FM <- na.omit(data[,c(var_subset, FM_test)])
    compare_clusters <- function(clus) {FM_index(clus, as.numeric(factor(Data_FM[,FM_test], levels=unique(Data_FM[,FM_test]))), assume_sorted_vectors = TRUE)}
    
    clusters_performance <- sapply(dend_clusters, compare_clusters)
    
    par(mar=c(4,4,2,2), mfrow=c(1,1))
    HC_clusters[["FM index testing"]] <- dotchart(sort(clusters_performance), xlim = FM_lim,
                                                  xlab = "Fowlkes-Mallows Index (from 0 to 1)",
                                                  main = "Performance of clustering algorithms",
                                                  pch = 19)
    
    
    if(export_plots=="png") {
      png(filename=paste(export_path, "FM Clustering Performance_", format(Sys.time(), "%Y-%m-%d %Hhr %Mmin %Ssec"),".png", sep=""),
          type="cairo", units="in",
          width=width, height=height,
          pointsize=point_size,
          res=dpi)
      
      par(mar=c(4,4,2,2), mfrow=c(1,1))
      HC_clusters[["FM index testing"]] <- dotchart(sort(clusters_performance), xlim = FM_lim,
                                                    xlab = "Fowlkes-Mallows Index (from 0 to 1)",
                                                    main = "Performance of clustering algorithms",
                                                    pch = 19)
      dev.off()
    }
  }
  
  if(export_plots=="png") {
    for (i in names(HC_clusters)) {
      png(filename=paste("HC_",i,"_", format(Sys.time(), "%Y-%m-%d %Hhr %Mmin %Ssec"), ".png", sep=""),
          type="cairo", units="in",
          width=width, height=height,
          pointsize=point_size,
          res=dpi)
      print(HC_clusters[[i]])
      dev.off()
    }
  } else if(export_plots=="pdf") {
    if(cluster_by=="variables") {
      pdf(paste0(export_path, "/HC_Clustering_by_Variables_", format(Sys.time(), "%Y-%m-%d %Hhr %Mmin %Ssec"),".pdf"),
          width=width, height=height)
    } else if(cluster_by=="observations") {
      pdf(paste0(export_path, "/HC_Clustering_by_Observations_", format(Sys.time(), "%Y-%m-%d %Hhr %Mmin %Ssec"),".pdf"),
          width=width, height=height)
    }
    
    print(HC_clusters)
    
    for(i in 1:length(Data_HC_clusters)) {
      dendextend_name <- paste("Base_",agglomeration_labels[i], sep="")
      HC_dendextend[[dendextend_name]] <- as.dendrogram(Data_HC_clusters[[i]])
      
      if(horizont==FALSE) {
        #pryr_plot1 %<a-% {
        par(mgp=c(1.6,0.4,0),mar=c(dend_mar[1],4,2,dend_mar[2]), mfrow=c(1,1))
        HC_dendextend[[dendextend_name]] %>%
          color_branches(.,k=k, col=cluster_colours) %>%
          set("labels_col", value=cluster_colours, k=k) %>%
          set("labels_cex", base_lab_size) %>%
          set("branches_lwd", base_lwd) %>% plot(horiz=horizont, yaxt="n", ylab=paste("Dissimilarity (", dist_measure,")", sep=""),
                                                 main=paste(dataset_label, " clusters (", agglomeration_labels[i], " linkage)", sep=""))
        
        axis(2, cex.axis=0.9, tck=-0.013) 
        
        if(draw_rect==TRUE) {
          HC_dendextend[[dendextend_name]] %>% rect.dendrogram(k=k, horiz = horizont,
                                                               border = "grey70", lty = 5, lwd = base_lwd)
        }
        
        HC_dendextend[[dendextend_name]] %<>% set("labels_col", value = cluster_colours, k= k)
        
        if(cluster_legend==TRUE) {
          # For the cluster legend...
          categories <- cutree(HC_dendextend[[dendextend_name]], k = k)
          colors_HC <- labels_colors(HC_dendextend[[dendextend_name]])[categories %>% sort %>% names]
          clust_labs <- colors_HC %>% unique
          clust_labs <- clust_labs[order(match(clust_labs, cluster_colours[1:length(clust_labs)]))]
          clust_leg <- unique(categories)
          
          if(!is.null(cluster_labels)) {
            clust_leg <- cluster_labels
          }
          legend("topright", legend = clust_leg, fill = clust_labs, cex = 0.9, bty="n")
        }
        #}
        #HC_clusters[[dendextend_name]] <- pryr_plot1
      } else if(horizont==TRUE) {
        #pryr_plot1 %<a-% {
        par(mgp=c(1.6,0.4,0),mar=c(dend_mar[1],2,2,dend_mar[2]), mfrow=c(1,1))
        HC_dendextend[[dendextend_name]] %>%
          color_branches(.,k=k, col=cluster_colours) %>%
          set("labels_col", value=cluster_colours, k=k) %>%
          set("labels_cex", base_lab_size) %>%
          set("branches_lwd", base_lwd) %>% plot(horiz=horizont, yaxt="n", xaxt="n", xlab=paste("Dissimilarity (", dist_measure,")", sep=""),
                                                 main=paste(dataset_label, " clusters (", agglomeration_labels[i], " linkage)", sep=""))
        
        axis(1, cex.axis=0.9, tck=-0.013) 
        
        if(draw_rect==TRUE) {
          HC_dendextend[[dendextend_name]] %>% rect.dendrogram(k=k, horiz = horizont,
                                                               border = "grey70", lty = 5, lwd = base_lwd)
        }
        
        HC_dendextend[[dendextend_name]] %<>% set("labels_col", value = cluster_colours, k = k)
        
        if(cluster_legend==TRUE) {
          # For the cluster legend...
          categories <- cutree(HC_dendextend[[dendextend_name]], k = k)
          colors_HC <- labels_colors(HC_dendextend[[dendextend_name]])[categories %>% sort %>% names]
          clust_labs <- colors_HC %>% unique
          clust_labs <- clust_labs[order(match(clust_labs, cluster_colours[1:length(clust_labs)]))]
          clust_leg <- unique(categories)
          
          if(!is.null(cluster_labels)) {
            clust_leg <- cluster_labels
          }
          legend("topleft", legend = clust_leg, fill = clust_labs, cex = 0.9, bty="n")
        }
      }
      #HC_clusters[[dendextend_name]] <- pryr_plot1
      #}
    }
    
    corrplot.mixed(HC_cor, upper="number", lower="pie", tl.col="black", tl.cex=0.8, number.cex=0.8, pch.cex=1.5, title=paste("Linkage Method Comparison (", dist_measure, " distance)", sep=""), mar=c(0,0,1,0))
    
    if(!is.null(FM_test) & cluster_by=="observations") {
      par(mar=c(4,4,2,2), mfrow=c(1,1))
      HC_clusters[["FM index testing"]] <- dotchart(sort(clusters_performance), xlim = FM_lim,
                                                    xlab = "Fowlkes-Mallows Index (from 0 to 1)",
                                                    main = "Performance of clustering algorithms",
                                                    pch = 19)
    }
    
    dev.off()
  }
  print(HC_clusters)
  return(list(Data_HC_clusters,HC_clusters))
}


# 3) Plot the clustering results against a temporal profile
Overlay_Clusters <- function(data, cluster_data, cluster_labels, var_subset, var_labels=var_subset, x_var, x_lab="X label", colours, facet_by="linkage method", cluster_type="arbitrary", 
                             k, export_plots="pdf", export_results=TRUE, export_path=getwd(), height=5, width=5, point_size=10, dpi=500) {
  
  #Install and load the required libraries 
  list.of.packages <- c("biotools", "asbio", "readr", "lattice", "corrplot", "fmsb", "latticeExtra", "RColorBrewer",
                        "reshape2", "gridExtra", "car", "psych", "data.table", "gridGraphics", "ggcorrplot", "corrgram",
                        "Hmisc", "xlsx", "factoextra", "NbClust", "dendextend", "fpc", "dtw", "dtwclust", "TSclust", 
                        "rlist", "pryr", "GGally", "ggpmisc", "lmtest", "cpm", "ecp", "changepoint", "magrittr", "scales")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  
  library(biotools)
  library(asbio)
  library(readr)
  library(lattice)
  library(corrplot)
  library(fmsb)
  library(latticeExtra)
  library(RColorBrewer)
  library(reshape2)
  library(gridExtra)
  library(car)
  library(psych)
  library(data.table)
  library(gridGraphics)
  library(ggcorrplot)
  library(corrgram)
  library(Hmisc)
  library(xlsx)
  library(factoextra)
  library(NbClust)
  library(dendextend)
  library(fpc)
  library(dtw)
  library(dtwclust)
  library(TSclust)
  library(rlist)
  library(pryr)
  library(GGally)
  library(ggpmisc)
  library(lmtest)
  library(cpm)
  library(ecp)
  library(changepoint)
  library(magrittr)
  library(scales)
  
  Data_List <- list()
  Plot_List <- list()
  
  Data <- na.omit(data)
  Data <- Data[,c(var_subset,x_var)]
  Data[,x_var] <- as.numeric(as.character(Data[,x_var]))
  
  if(cluster_type=="HC" & is.list(cluster_data)) {
    for(i in names(cluster_data)) {
      clustname <- paste(i, sep="")
      clust_data <- cutree(as.dendrogram(cluster_data[[i]]), k=k)
      Data_List[[i]] <- clust_data
    }
    plot_lnm <- names(Data_List)
    
  } else if(cluster_type!="HC") {
    if(is.list(cluster_data)) { for(i in seq_along(cluster_data)) { Data_List[[i]] <- cluster_data } } else Data_List <- cluster_data
    if(is.list(cluster_data)) plot_lnm <- names(cluster_data) else plot_lnm <- colnames(Data_List)
  }
  
  Plot_Data <- cbind(Data, as.data.frame(Data_List))
  
  if(export_results==TRUE) {
    fwrite(Plot_Data, paste0(export_path, "/HC_clusters_", format(Sys.time(), "%Y-%m-%d %Hhr %Mmin %Ssec"), ".csv"), na = NA)
  }
  
  if(facet_by=="none") {
    for(i in var_subset) {
      for(j in plot_lnm) {
        pltName <- paste(i,"_",j," Linkage", sep="")
        #graph_data[[i]]["Class"] <- factor(graph_data[[i]]["Class"], levels=unique(graph_data[[1]]))
        Plot_List[[pltName]] <- ggplot(data=Plot_Data, aes_string(x=x_var, y=c(i), colour=factor(Plot_Data[,j]), fill=factor(Plot_Data[,j]), shape=factor(Plot_Data[,j]))) +
          #geom_ribbon(data=graph_data[[i]], aes(x=Age/1000, ymin=min_err,ymax=max_err), colour="grey50", fill=c("#00AFBB"), alpha=0.1, inherit.aes=FALSE) +
          #geom_hline(colour="grey55", lty=2, lwd=0.7, yintercept = c(10,50)) +
          geom_path(colour="grey15", lwd=0.7) + # geom_path follows the order of data passed to it
          geom_point(colour="black", size=2.5) +
          #geom_smooth(method="lm", formula=my.formula,show.legend=FALSE, size=0.4) +
          scale_shape_manual(values=c(21:25, 0:20), name=NULL, labels=cluster_labels) +
          scale_colour_manual(values=colours) +
          scale_fill_manual(values=colours, name=NULL, labels=cluster_labels) +
          scale_x_continuous(breaks=pretty_breaks(n=10)) +
          #scale_y_continuous(breaks=pretty_breaks(n=10), limits=c(NA,100)) +
          theme(aspect.ratio = 1,
                axis.text.x = element_text(colour="black", angle = 0, size=10, margin=margin(c(0.1,0.1,0.1,0.1), unit="cm")),
                axis.text.y = element_text(colour="black", size=10, margin=margin(c(0.1,0.1,0.1,0.1), unit="cm")),
                axis.ticks.length=unit(0.1, "cm"),
                #axis.ticks=element_line(),
                plot.margin = grid::unit(c(1, 1, 1, 1), "lines"),
                axis.line=element_line(colour="grey75"),
                panel.grid.minor=element_blank(),
                panel.grid.major.x=element_blank(),
                #panel.grid.major.y=element_line(colour="grey65",linetype="solid"),
                panel.grid.major=element_blank(),
                panel.border=element_rect(fill=NA, colour="grey40"),
                panel.background=element_blank(),
                axis.title.x=element_text(colour="black", margin=margin(t=0.1, unit="cm")),
                axis.title.y=element_text(colour="black", margin=margin(r=0.1, unit="cm")),
                legend.position="bottom", legend.direction="horizontal",
                legend.key=element_blank(),
                legend.margin=margin(t=0, unit="cm"),
                legend.text=element_text(size=11),
                legend.text.align=0,
                legend.title=element_text(face="bold")) + 
          labs(x=x_lab, y=var_labels[match(i, var_subset)])#expression(P[III]*IP[25]*"-derived SpSIC (%)"), fill=NULL)
      }
    }
  } else if(facet_by=="linkage method") {
    Plot_Data <- melt(setDT(Plot_Data), id.vars=c(x_var, plot_lnm), measure.vars=c(var_subset))
    Plot_Data <- as.data.frame(melt(Plot_Data, id.vars=c(x_var, "value", "variable")))
    Plot_Data[,"value.1"] <- factor(Plot_Data[,"value.1"], levels=unique(Plot_Data[,"value.1"]))
    for(i in c(var_subset)) {
      pltName <- paste(c(i),"_Linkage Method faceting", sep="")
      Plot_List[[pltName]] <- ggplot(data=Plot_Data[Plot_Data[,c("variable")] %in% c(i),], aes_string(x=x_var, y="value", colour="value.1", fill="value.1", shape="value.1")) +
        #geom_ribbon(data=graph_data[[i]], aes(x=Age/1000, ymin=min_err,ymax=max_err), colour="grey50", fill=c("#00AFBB"), alpha=0.1, inherit.aes=FALSE) +
        #geom_hline(colour="grey55", lty=2, lwd=0.7, yintercept = c(10,50)) +
        geom_path(colour="grey15", group=1, lwd=0.7) + # geom_path follows the order of data passed to it
        geom_point(colour="black", size=2.5) +
        #geom_smooth(method="lm", formula=my.formula,show.legend=FALSE, size=0.4) +
        facet_wrap(~variable.1, scales="free_y") +
        scale_shape_manual(values=c(21:25, 0:20), name=NULL, labels=cluster_labels) +
        scale_colour_manual(values=colours) +
        scale_fill_manual(values=colours, name=NULL, labels=cluster_labels) +
        scale_x_continuous(breaks=pretty_breaks(n=10)) +
        #scale_y_continuous(breaks=pretty_breaks(n=10), limits=c(NA,100)) +
        theme(aspect.ratio = 1,
              axis.text.x = element_text(colour="black", angle = 0, size=10, margin=margin(c(0.1,0.1,0.1,0.1), unit="cm")),
              axis.text.y = element_text(colour="black", size=10, margin=margin(c(0.1,0.1,0.1,0.1), unit="cm")),
              axis.ticks.length=unit(0.1, "cm"),
              #axis.ticks=element_line(),
              plot.margin = grid::unit(c(1, 1, 1, 1), "lines"),
              axis.line=element_line(colour="grey75"),
              panel.grid.minor=element_blank(),
              panel.grid.major.x=element_blank(),
              #panel.grid.major.y=element_line(colour="grey65",linetype="solid"),
              panel.grid.major=element_blank(),
              panel.border=element_rect(fill=NA, colour="grey40"),
              panel.background=element_blank(),
              axis.title.x=element_text(colour="black", margin=margin(t=0.1, unit="cm")),
              axis.title.y=element_text(colour="black", margin=margin(r=0.1, unit="cm")),
              legend.position="bottom", legend.direction="horizontal",
              legend.key=element_blank(),
              legend.margin=margin(t=0, unit="cm"),
              legend.text=element_text(size=11),
              legend.text.align=0,
              legend.title=element_text(face="bold")) + 
        labs(x=x_lab, y=var_labels[match(c(i), var_subset)], colour="Clusters")#expression(P[III]*IP[25]*"-derived SpSIC (%)"), fill=NULL)
    }
  } else if(facet_by=="variables") {
    Plot_Data <- melt(setDT(Plot_Data), id.vars=c(x_var, plot_lnm), measure.vars=c(var_subset))
    Plot_Data <- as.data.frame(melt(Plot_Data, id.vars=c(x_var, "value", "variable")))
    Plot_Data[,"value.1"] <- factor(Plot_Data[,"value.1"], levels=unique(Plot_Data[,"value.1"]))
    for(i in plot_lnm) {
      pltName <- paste(c(i)," Linkage_All Variables", sep="")
      Plot_List[[pltName]] <- ggplot(data=Plot_Data[Plot_Data[,c("variable.1")] %in% c(i),], aes_string(x=x_var, y="value", colour="value.1", fill="value.1", shape="value.1")) +
        #geom_ribbon(data=graph_data[[i]], aes(x=Age/1000, ymin=min_err,ymax=max_err), colour="grey50", fill=c("#00AFBB"), alpha=0.1, inherit.aes=FALSE) +
        #geom_hline(colour="grey55", lty=2, lwd=0.7, yintercept = c(10,50)) +
        geom_path(colour="grey15", group=1, lwd=0.7) + # geom_path follows the order of data passed to it
        geom_point(colour="black", size=2.5) +
        #geom_smooth(method="lm", formula=my.formula,show.legend=FALSE, size=0.4) +
        facet_wrap(~variable, scales="free_y") +
        scale_shape_manual(values=c(21:25, 0:20), name=NULL, labels=cluster_labels) +
        scale_colour_manual(values=colours) +
        scale_fill_manual(values=colours, name=NULL, labels=cluster_labels) +
        scale_x_continuous(breaks=pretty_breaks(n=10)) +
        #scale_y_continuous(breaks=pretty_breaks(n=10), limits=c(NA,100)) +
        theme(aspect.ratio=1,
              axis.text.x = element_text(colour="black", angle = 0, size=10, margin=margin(c(0.1,0.1,0.1,0.1), unit="cm")),
              axis.text.y = element_text(colour="black", size=10, margin=margin(c(0.1,0.1,0.1,0.1), unit="cm")),
              axis.ticks.length=unit(0.1, "cm"),
              #axis.ticks=element_line(),
              plot.margin = grid::unit(c(1, 1, 1, 1), "lines"),
              axis.line=element_line(colour="grey75"),
              panel.grid.minor=element_blank(),
              panel.grid.major.x=element_blank(),
              #panel.grid.major.y=element_line(colour="grey65",linetype="solid"),
              panel.grid.major=element_blank(),
              panel.border=element_rect(fill=NA, colour="grey40"),
              panel.background=element_blank(),
              axis.title.x=element_text(colour="black", margin=margin(t=0.1, unit="cm")),
              axis.title.y=element_text(colour="black", margin=margin(r=0.1, unit="cm")),
              legend.position="bottom", legend.direction="horizontal",
              legend.key=element_blank(),
              legend.margin=margin(t=0, unit="cm"),
              legend.text=element_text(size=11),
              legend.text.align=0,
              legend.title=element_text(face="bold")) + 
        labs(x=x_lab, y=c(i), colour="Clusters")#expression(P[III]*IP[25]*"-derived SpSIC (%)"), fill=NULL)
    }
  }
  
  if(export_plots=="png") {
    for (i in names(Plot_List)) {
      png(filename=paste0(export_path, "/", i, "_", format(Sys.time(), "%Y-%m-%d %Hhr %Mmin %Ssec"), ".png"),
          type="cairo", units="in",
          width=width, height=height,
          pointsize=point_size,
          res=dpi)
      print(Plot_List[[i]])
      dev.off()
    }
  } else if(export_plots=="pdf") {
    if(facet_by=="linkage method") {
      pdf(paste0(export_path, "/Cluster_Overlay_by Linkage Method_", format(Sys.time(), "%Y-%m-%d %Hhr %Mmin %Ssec"), ".pdf"),
          width=width, height=height)
    } else if(facet_by=="variables") {
      pdf(paste0(export_path, "/Cluster_Overlay_by Variable_", format(Sys.time(), "%Y-%m-%d %Hhr %Mmin %Ssec"), ".pdf"),
          width=width, height=height)
    } else if(facet_by=="none") {
      pdf(paste0(export_path, "/Cluster_Overlay_Separate Plots_", format(Sys.time(), "%Y-%m-%d %Hhr %Mmin %Ssec"), ".pdf"),
          width=width, height=height)
    }
    print(Plot_List)
    dev.off()
  }
  print(Plot_List)
  return(Plot_List)
}

# 4) Align two time series from e.g. downcore records using DTW techniques
DTW_alignment <- function(query_data, ref_data, query_title="Query", ref_title="Reference", query_subset=NULL, ref_subset=NULL, z_norm=list(TRUE, FALSE, FALSE), reinterp=FALSE,
                          x_var, step_pattern="symmetric2", window_type="sakoechiba", window_size=20, MVM_elasticity=50, RJ_step_settings=c(4, "d", FALSE), open_begin=FALSE, open_end=FALSE, x_align=TRUE,
                          y_offset=NULL, x_label=x_var, y_labels_query=waiver(), y_labels_ref=waiver(), colours=brewer.pal(9, "Set1"), match_subset=NULL, match_min=NULL, sample_freq=1, match_vis=0.8, grang_order=3,
                          x_rounding_order=-1, y_rounding_order=-1, export_plots="pdf", width=10, height=5, point_size=10,
                          dpi=500, export_results=TRUE, export_path=getwd(), dtw_distance="Euclidean") {
  
  #Install and load the required libraries 
  list.of.packages <- c("biotools", "asbio", "readr", "lattice", "corrplot", "fmsb", "latticeExtra", "RColorBrewer",
                        "reshape2", "gridExtra", "car", "psych", "data.table", "gridGraphics", "ggcorrplot", "corrgram",
                        "Hmisc", "xlsx", "factoextra", "NbClust", "dendextend", "fpc", "dtw", "dtwclust", "TSclust", 
                        "rlist", "pryr", "GGally", "ggpmisc", "lmtest", "cpm", "ecp", "changepoint", "magrittr", "scales", "BBmisc")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  
  library(biotools)
  library(asbio)
  library(readr)
  library(lattice)
  library(corrplot)
  library(fmsb)
  library(latticeExtra)
  library(RColorBrewer)
  library(reshape2)
  library(gridExtra)
  library(car)
  library(psych)
  library(data.table)
  library(gridGraphics)
  library(ggcorrplot)
  library(corrgram)
  library(Hmisc)
  library(xlsx)
  library(factoextra)
  library(NbClust)
  library(dendextend)
  library(fpc)
  library(dtw)
  library(dtwclust)
  library(TSclust)
  library(rlist)
  library(pryr)
  library(GGally)
  library(ggpmisc)
  library(lmtest)
  library(cpm)
  library(ecp)
  library(changepoint)
  library(magrittr)
  library(scales)
  library(BBmisc)
  
  #Modified plotting functions from the DTW package
  Threeway_Alt <- function (d, xts = NULL, yts = NULL, time_x = NULL, time_y = NULL, type.align = "l", type.ts = "l", 
                            match.indices = NULL, top_match = NULL, margin_def = 4, inner.margin = 0.2, title.margin = 1.5, 
                            xlab = "Query index", ylab = "Reference index", main = "Timeseries alignment",
                            colours = brewer.pal(9, "Set1"), ...) 
  {
    if (is.null(xts) || is.null(yts)) {
      xts <- d$query
      yts <- d$reference
    }
    if (is.null(xts) || is.null(yts)) 
      stop("Original timeseries are required")
    xts <- as.matrix(xts)
    yts <- as.matrix(yts)
    
    if(is.null(time_x) || is.null(time_y)) {
      stop("Time series time_x and time_y must be provided!")
      #print("Custom x-axes for time series not provided! Using equi-spaced index...")
      #time_x <- d$index1
      #time_y <- d$index2
    }
    
    if (ncol(xts) > 1 || ncol(yts) > 1) 
      stop("Only single-variate timeseries can be displayed. (You may want to extract a column for visualization purposes.)")
    def.par <- par(no.readonly = TRUE)
    layout(matrix(c(3, 1, 0, 2), 2, 2, byrow = TRUE), c(1, 3), 
           c(3, 1), TRUE)
    imar <- inner.margin
    bmar <- margin_def
    lmar <- margin_def
    tmar <- margin_def + title.margin
    rmar <- margin_def
    mlab = margin_def/2
    mtex = margin_def/6
    nn <- length(xts)
    mm <- length(yts)
    par(mar = c(imar, imar, tmar, rmar))
    plot(d$index1, d$index2, type = type.align, xlim = c(1, nn), 
         ylim = c(1, mm), ax = FALSE, main = main, ...)
    if (length(match.indices) == 1 || is.null(match.indices)) {
      match.indices <- seq(from = 1, to = length(d$index1), 
                           length.out = length(d$index1))
    }
    if (!is.null(match.indices)) {
      idx <- match.indices
      
      if(is.null(idx)) idx <- seq(from = 1, to = length(d$index1),
                                  length.out = match.indices)
      
      segments(d$index1[idx], 0, d$index1[idx], d$index2[idx], 
               col = colours[3], lty = 3)
      segments(0, d$index2[idx], d$index1[idx], d$index2[idx], 
               col = colours[3], lty = 3)
    }
    
    if(!is.null(top_match)) {
      top_idx <- as.numeric(top_match)
      segments(d$index1[top_idx], 0, d$index1[top_idx], d$index2[top_idx], 
               col = colours[4], lty = 5, lwd=1.3)
      segments(0, d$index2[top_idx], d$index1[top_idx], d$index2[top_idx], 
               col = colours[4], lty = 5, lwd=1.3)
    }
    
    box()
    par(mar = c(bmar, imar, imar, rmar))
    plot(xts ~ time_x, type = type.ts, xlab = xlab, mgp = c(mlab, 
                                                            mtex, 0), ax = FALSE, col = colours[1])
    axis(1)
    axis(2)
    box()
    par(mar = c(imar, lmar, tmar, imar))
    plot(time_y ~ yts, xlim = rev(range(yts[!is.na(yts)])), type = type.ts, 
         ylab = ylab, mgp = c(mlab, mtex, 0), ax = FALSE, col = colours[2])
    axis(3)
    axis(2)
    box()
    par(def.par)
  }
  
  #Lists to hold DTW results
  alignment_list <- list()
  DTW_Data_List <- list()
  DTW_Query_List <- list()
  DTW_Ref_List <- list()
  DTW_DistRank <- list()
  
  granger_results <- as.data.frame(matrix(c("DF difference","Wald Statistic","P-value", "Cumulative DTW distance", "Normalized DTW distance")))
  
  if(length(query_subset)!=length(ref_subset)) {
    stop("The lengths of query_subset and ref_subset must be equal! Please adjust.")
  }
  
  if(is.character(query_subset) & is.character(ref_subset)) {
    que_data <- query_data[c(x_var, query_subset)]
  } else if(!is.character(query_subset) & is.character(ref_subset)) {
    cat("Argument query_subset not provided! Assuming equality to ref_subset...")
    query_subset <- ref_subset
  } else if(is.character(query_subset) & !is.character(ref_subset)) {
    cat("Argument ref_subset not provided! Assuming equality to query_subset...")
    ref_subset <- query_subset
  } else if(!is.character(ref_subset) & !is.character(query_subset)) {
    stop("The query_subset and/or ref_subset arguments must be character values! Please adjust.")
  }
  
  for(i in seq_along(query_subset)) {
    if(reinterp==FALSE) {
      xdata <- que_data[[c(x_var)]]
      qu <- que_data[[query_subset[i]]]
      ref <- ref_data[[ref_subset[i]]]
    } else if(reinterp==TRUE) {
      xdata <- reinterpolate(que_data[[c(x_var)]], new.length=length(ref_data[[c(x_var)]]))
      qu <- reinterpolate(que_data[[query_subset[i]]], new.length=length(ref_data[[ref_subset[i]]]))
      ref <- ref_data[[ref_subset[i]]]
    }
    
    if(isTRUE(z_norm[[1]])) {
      qu <- scale(qu)
      ref <- scale(ref)
    }
    
    DTWname <- paste("DTW alignment using ", query_subset[i], query_title, " and ", ref_subset[i], " ", ref_title, sep="") #rabinerJuangStepPattern(type,slope.weighting,smoothed)
    if(step_pattern=="RJ") {
      alignment_list[[DTWname]] <- try(dtw(qu, ref, dist.method=dtw_distance, step.pattern = rabinerJuangStepPattern(as.integer(RJ_step_settings[1]), as.character(RJ_step_settings[2]), RJ_step_settings[3]), open.end=open_end, open.begin=open_begin, window.type=window_type, window.size=window_size, k=TRUE))
    } else if(step_pattern=="MVM") {
      alignment_list[[DTWname]] <- try(dtw(qu, ref, dist.method=dtw_distance, step.pattern = mvmStepPattern(elasticity=MVM_elasticity), open.end=open_end, open.begin=open_begin, window.type=window_type, window.size=window_size, k=TRUE))
    } else if(step_pattern!="MVM" & step_pattern!="RJ") {
      alignment_list[[DTWname]] <- try(dtw(qu, ref, dist.method=dtw_distance, step.pattern = eval(parse(text=step_pattern)), open.end=open_end, open.begin=open_begin, window.type=window_type, window.size=window_size, k=TRUE))
    }
    if(inherits(alignment_list[[DTWname]], "try-error")) {
      stop("Optimal warping path not found under the given constraints! Try changing step_pattern, open_end, open_begin.")
    }
    #dtwPlot(alignment_list[[DTWname]], type="twoway", off=y_offset, match.indices = matching_freq)
    if(!isTRUE(z_norm[[1]])) {
      DTW_Data_List[[paste(query_subset[i], "vs", ref_subset[i])]] <- cbind(index1=c(alignment_list[[DTWname]]$index1), index2=c(alignment_list[[DTWname]]$index2), Depth1=as.data.frame(xdata[as.numeric(alignment_list[[DTWname]]$index1)]), Depth2=as.data.frame(ref_data[[x_var]][as.numeric(alignment_list[[DTWname]]$index2)]),
                                                                            var1=as.data.frame(qu[as.numeric(alignment_list[[DTWname]]$index1)]), var2=as.data.frame(ref[as.numeric(alignment_list[[DTWname]]$index2)]))
    } else if(isTRUE(z_norm[[1]])) {
      if(isTRUE(reinterp)) {
        quote_norm <- reinterpolate(que_data[[query_subset[i]]], new.length=length(ref_data[[ref_subset[i]]]))
      } else if(!isTRUE(reinterp)) {
        quote_norm <- que_data[[query_subset[i]]]
      }
      ref_norm <- ref_data[[ref_subset[i]]]
      DTW_Data_List[[paste(query_subset[i], "vs", ref_subset[i])]] <- cbind(index1=c(alignment_list[[DTWname]]$index1), index2=c(alignment_list[[DTWname]]$index2), Depth1=as.data.frame(xdata[as.numeric(alignment_list[[DTWname]]$index1)]), Depth2=as.data.frame(ref_data[[x_var]][as.numeric(alignment_list[[DTWname]]$index2)]),
                                                                            var1=as.data.frame(quote_norm[as.numeric(alignment_list[[DTWname]]$index1)]), var2=as.data.frame(ref_norm[as.numeric(alignment_list[[DTWname]]$index2)]))  
    }
    colnames(DTW_Data_List[[paste(query_subset[i], "vs", ref_subset[i])]]) <- c(paste0(query_title, "_index_query"), paste0(ref_title, "_index_ref"), paste(x_var, query_title, "query", sep="_"), paste(x_var, ref_title, "ref", sep="_"), paste(query_subset[i], query_title, "query", sep="_"), paste(ref_subset[i], ref_title, "ref", sep="_"))
    
    #Retrieving local costs (NOT cumulative) from DTW object (EXAMPLE - implement into DTW function!)
    DTW_DistRes <- alignment_list[[DTWname]]
    index_qu <- unlist(DTW_DistRes["index1"])
    index_ref <- unlist(DTW_DistRes["index2"])
    
    if(!identical(length(index_qu), length(index_ref))) {
      stop("The lengths of DTW indices are unequal! Aborting...")
    }
    
    local_costvec <- sapply(1:length(index_qu), function(x) DTW_DistRes[["localCostMatrix"]][index_qu[x],index_ref[x]])
    cumul_costvec <- sapply(1:length(index_qu), function(x) DTW_DistRes[["costMatrix"]][index_qu[x],index_ref[x]])
    DTW_Data_List[[paste(query_subset[i], "vs", ref_subset[i])]] <- cbind(DTW_Data_List[[paste(query_subset[i], "vs", ref_subset[i])]], Local_Cost=local_costvec, Cumul_Cost=cumul_costvec, Steps=c(DTW_DistRes[["stepsTaken"]], 0))
    
    names(local_costvec) <- names(cumul_costvec) <- 1:length(local_costvec)
    if(!is.null(match_min) & is.numeric(match_min)) {
      local_maxind <- as.numeric(names(local_costvec)[order(local_costvec, decreasing = FALSE)][1:match_min])
      DTW_DistRank[[i]] <- cbind(DTW_Data_List[[paste(query_subset[i], "vs", ref_subset[i])]][local_maxind,], Dist_RowIndex=local_maxind) 
    }
    
    #Reference and query data lists (for plots)
    DTW_Query_List[[paste(query_subset[i], "vs", ref_subset[i])]] <- cbind(as.data.frame(que_data)[,c(x_var, query_subset[i])])
    DTW_Ref_List[[paste(query_subset[i], "vs", ref_subset[i])]] <- cbind(as.data.frame(ref_data)[,c(x_var, ref_subset[i])])
    colnames(DTW_Query_List[[paste(query_subset[i], "vs", ref_subset[i])]]) <- c(x_var, query_subset[i])
    colnames(DTW_Ref_List[[paste(query_subset[i], "vs", ref_subset[i])]]) <- c(x_var, ref_subset[i])
    
    cumul_dist <- alignment_list[[DTWname]]$distance
    norm_dist <- alignment_list[[DTWname]]$normalizedDistance
    grang <- suppressWarnings(try(grangertest(que_data[[query_subset[i]]], ref_data[[ref_subset[i]]], order=grang_order)))
    
    if(!inherits(grang, "try-error")) {
      granger_results <- cbind(granger_results, c(grang$Df[2],grang$`F`[2], grang$`Pr(>F)`[2], cumul_dist, norm_dist))
    } else if(inherits(grang, "try-error")) {
      granger_results <- cbind(c("Cumulative Distance", "Normalised Distance"), cumul_dist, norm_dist)
    }
    
  }
  #return(list(DTW_Data_List, DTW_Query_List, DTW_Ref_List))}
  #colnames(granger_results) <- c("Parameter", paste(names(DTW_Data_List)))
  
  if(export_results==TRUE) {
    DTW_Results <- createWorkbook()
    for(i in seq_along(DTW_Data_List)) {
      sheetname <- paste(names(DTW_Data_List[i]),"DTW alignment")
      sheet <- createSheet(DTW_Results, paste(sheetname))
      addDataFrame(DTW_Data_List[[i]], sheet=sheet, startColumn=1, row.names=FALSE)
      
      if(!is.null(match_min) & is.numeric(match_min)) {
        addDataFrame(DTW_DistRank[[i]], sheet=sheet, startColumn=ncol(DTW_Data_List[[i]])+2, row.names=FALSE)
      }
    }
    #sheet_granger <- createSheet(DTW_Results, "Granger Test and DTW distances")
    #addDataFrame(granger_results, sheet=sheet_granger, row.names=FALSE)
    saveWorkbook(DTW_Results,paste0(export_path,"/DTW Alignment Results_", format(Sys.time(), "%Y-%m-%d %Hhr %Mmin %Ssec"), ".xlsx"))
  }
  
  GG_Plot_List <- list()
  
  #Boundary between data processing and plotting here!
  #return(DTW_Data_List)
  
  for(i in seq_along(DTW_Data_List)) {
    plot_segment <- DTW_Data_List[[i]]
    plot_query <- DTW_Query_List[[i]]
    plot_ref <- DTW_Ref_List[[i]]
    
    if(!is.numeric(z_norm[[2]]) & !is.numeric(z_norm[[3]])) {
      if(y_offset[1]!=0 & !is.null(y_offset[1])) {
        y_offset_value <- as.numeric(y_offset[1])*mean(plot_segment[,5])
        plot_segment[,6] <- plot_segment[,6] + y_offset_value
        plot_ref[,2] <- plot_ref[,2] + y_offset_value
        
        if(!is.null(match_min) & is.numeric(match_min)) {
          #y_distrank <-y_offset*mean(DTW_DistRank[[i]][,5])
          DTW_DistRank[[i]][,6] <- DTW_DistRank[[i]][,6] + y_offset_value
        }
      }
    } else if(is.numeric(z_norm[[2]]) & is.numeric(z_norm[[3]])) {
      y_offset <- c(0, "none")
      
      library(BBmisc)
      
      plot_query[,2] <- BBmisc::normalize(plot_query[,2], method="range", range=z_norm[[2]])
      plot_ref[,2] <- BBmisc::normalize(plot_ref[,2], method="range", range=z_norm[[3]])
      plot_segment[,5] <- BBmisc::normalize(plot_segment[,5], method="range", range=z_norm[[2]])
      plot_segment[,6] <- BBmisc::normalize(plot_segment[,6], method="range", range=z_norm[[3]])
      
      if(!is.null(match_min) & is.numeric(match_min)) {
        DTW_DistRank[[i]][,5] <- plot_segment[as.numeric(DTW_DistRank[[i]][,10]),5]
        DTW_DistRank[[i]][,6] <- plot_segment[as.numeric(DTW_DistRank[[i]][,10]),6]
      }
    }
    
    
    if(x_align==TRUE) {
      x_align_factor <- max(plot_ref[,x_var])/max(plot_query[,x_var])
      plot_segment[,3] <- plot_segment[,3] * x_align_factor
      plot_query[,x_var] <- plot_query[,x_var] * x_align_factor
      
      if(!is.null(match_min) & is.numeric(match_min)) {
        DTW_DistRank[[i]][,3] <- DTW_DistRank[[i]][,3] * x_align_factor
      }
    }
    
    if(!is.null(match_subset)) {
      if(x_align==TRUE) {
        match_sub <- match_subset*x_align_factor
      } else if(!isTRUE(x_align)) {
        match_sub <- as.numeric(match_subset)
      }
      plot_seg <- plot_segment[plot_segment[,3] %in% match_sub,]
    } else if(is.null(match_subset)) {
      plot_seg <- plot_segment
    }
    
    if(is.null(y_offset[1]) | y_offset[1]==0| y_offset[2]!="new_breaks") {
      plotname <- paste(names(DTW_Data_List)[i], " DTW Plot", sep="")
      GG_Plot_List[[plotname]] <- ggplot() +
        geom_line(data=plot_query, aes_string(x=c(x_var), y=colnames(plot_query)[2]), colour=colours[1]) +
        geom_line(data=plot_ref, aes_string(x=c(x_var), y=colnames(plot_ref)[2]), colour=colours[2]) +
        geom_segment(data=plot_seg[seq(1, nrow(plot_seg), sample_freq),], aes_string(x=colnames(plot_seg)[3], y=colnames(plot_seg)[5], xend=colnames(plot_seg)[4], yend=colnames(plot_seg)[6]), colour=colours[3], alpha=match_vis, linetype=3) +
        scale_y_continuous(name = paste(y_labels_query[i], "for", query_title, "(Query)"), 
                           sec.axis = sec_axis(~., name = paste(y_labels_ref[i], "for", query_title, "(Reference)"))) +
        theme(#aspect.ratio = 1,
          axis.text.x = element_text(colour="black", angle = 0, size=10, margin=margin(c(0.1,0.1,0.1,0.1), unit="cm")),
          axis.text.y = element_text(colour="black", size=10, margin=margin(c(0.1,0.1,0.1,0.1), unit="cm")),
          axis.ticks.length=unit(0.1, "cm"),
          #axis.ticks=element_line(),
          plot.margin = grid::unit(c(1, 1, 1, 1), "lines"),
          axis.line=element_line(colour="grey75"),
          panel.grid.minor=element_blank(),
          panel.grid.major.x=element_blank(),
          #panel.grid.major.y=element_line(colour="grey65",linetype="solid"),
          panel.grid.major=element_blank(),
          panel.border=element_rect(fill=NA, colour="grey40"),
          panel.background=element_blank(),
          axis.title.x=element_text(colour="black", margin=margin(t=0.1, unit="cm")),
          axis.title.y=element_text(colour="black", margin=margin(r=0.1, unit="cm")),
          legend.position="bottom", legend.direction="horizontal",
          legend.key=element_blank(),
          legend.margin=margin(t=0, unit="cm"),
          legend.text=element_text(size=11),
          legend.text.align=0,
          legend.title=element_text(face="bold")) +
        labs(x=x_label, y=paste(y_labels_query[i], "for", query_title, "(Query)"))
    } else if(y_offset[1]!=0 & !is.null(y_offset[1]) & !isTRUE(x_align) & y_offset[2]=="new_breaks") { #| y_offset!=0 & isTRUE(x_align) & isTRUE(reinterp)
      plotname <- paste(names(DTW_Data_List)[i], " DTW Plot", sep="")
      GG_Plot_List[[plotname]] <- ggplot() +
        geom_line(data=plot_query, aes_string(x=c(x_var), y=colnames(plot_query)[2]), colour=colours[1]) +
        geom_line(data=plot_ref, aes_string(x=c(x_var), y=colnames(plot_ref)[2]), colour=colours[2]) +
        geom_segment(data=plot_seg[seq(1, nrow(plot_seg), sample_freq),], aes_string(x=colnames(plot_seg)[3], y=colnames(plot_seg)[5], xend=colnames(plot_seg)[4], yend=colnames(plot_seg)[6]), colour=colours[3], alpha=match_vis, linetype=3) +
        scale_y_continuous(name = paste(y_labels_query[i], "for", query_title, "(Query)"), 
                           sec.axis = sec_axis(~ . - y_offset_value, name = paste(y_labels_ref[i], "for", ref_title, "(Reference)"), 
                                               breaks = seq(0,round(max(plot_seg*11.5/10),y_rounding_order),round((max(plot_seg*11.5/10)/5),y_rounding_order)))) +
        theme(#aspect.ratio = 1,
          axis.text.x = element_text(colour="black", angle = 0, size=10, margin=margin(c(0.1,0.1,0.1,0.1), unit="cm")),
          axis.text.y = element_text(colour="black", size=10, margin=margin(c(0.1,0.1,0.1,0.1), unit="cm")),
          axis.ticks.length=unit(0.1, "cm"),
          #axis.ticks=element_line(),
          plot.margin = grid::unit(c(1, 1, 1, 1), "lines"),
          axis.line=element_line(colour="grey75"),
          panel.grid.minor=element_blank(),
          panel.grid.major.x=element_blank(),
          #panel.grid.major.y=element_line(colour="grey65",linetype="solid"),
          panel.grid.major=element_blank(),
          panel.border=element_rect(fill=NA, colour="grey40"),
          panel.background=element_blank(),
          axis.title.x=element_text(colour="black", margin=margin(t=0.1, unit="cm")),
          axis.title.y=element_text(colour="black", margin=margin(r=0.1, unit="cm")),
          legend.position="bottom", legend.direction="horizontal",
          legend.key=element_blank(),
          legend.margin=margin(t=0, unit="cm"),
          legend.text=element_text(size=11),
          legend.text.align=0,
          legend.title=element_text(face="bold")) +
        labs(x=x_label, y=paste(y_labels_query[i], "for", query_title, "(Query)"))
    } else if(y_offset[1]!=0 & y_offset[2]=="new_breaks" & isTRUE(x_align)) { #& !isTRUE(reinterp)
      plotname <- paste(names(DTW_Data_List)[i], " DTW Plot", sep="")
      GG_Plot_List[[plotname]] <- ggplot() +
        geom_line(data=plot_query, aes_string(x=c(x_var), y=colnames(plot_query)[2]), colour=colours[1]) +
        geom_line(data=plot_ref, aes_string(x=c(x_var), y=colnames(plot_ref)[2]), colour=colours[2]) +
        geom_segment(data=plot_seg[seq(1, nrow(plot_seg), sample_freq),], aes_string(x=colnames(plot_seg)[3], y=colnames(plot_seg)[5], xend=colnames(plot_seg)[4], yend=colnames(plot_seg)[6]), colour=colours[3], alpha=match_vis, linetype=3) +
        scale_y_continuous(name = paste(y_labels_query[i], "for",query_title, "(Query)"), 
                           sec.axis = sec_axis(~ . - y_offset_value, name = paste(y_labels_ref[i], "for", ref_title, "(Reference)"), 
                                               breaks = seq(0,round(max(plot_seg[,6]*11.5/10),y_rounding_order),round((max(plot_seg[,6]*11.5/10)/5),y_rounding_order)))) +
        scale_x_continuous(position="top", name = paste(x_label, "for", ref_title, "(Reference)"), 
                           sec.axis = sec_axis(~ ./x_align_factor, name = paste(x_label, "for", query_title, "(Query)"), 
                                               breaks = seq(0,round(max(plot_seg[,3]*11.5/10),x_rounding_order),round((max(plot_seg[,3]*11.5/10)/5),x_rounding_order)))) +
        theme(#aspect.ratio = 1,
          axis.text.x = element_text(colour="black", angle = 0, size=10, margin=margin(c(0.1,0.1,0.1,0.1), unit="cm")),
          axis.text.y = element_text(colour="black", size=10, margin=margin(c(0.1,0.1,0.1,0.1), unit="cm")),
          axis.ticks.length=unit(0.1, "cm"),
          #axis.ticks=element_line(),
          plot.margin = grid::unit(c(1, 1, 1, 1), "lines"),
          axis.line=element_line(colour="grey75"),
          panel.grid.minor=element_blank(),
          panel.grid.major.x=element_blank(),
          #panel.grid.major.y=element_line(colour="grey65",linetype="solid"),
          panel.grid.major=element_blank(),
          panel.border=element_rect(fill=NA, colour="grey40"),
          panel.background=element_blank(),
          axis.title.x=element_text(colour="black", margin=margin(t=0.1, unit="cm")),
          axis.title.y=element_text(colour="black", margin=margin(r=0.1, unit="cm")),
          legend.position="bottom", legend.direction="horizontal",
          legend.key=element_blank(),
          legend.margin=margin(t=0, unit="cm"),
          legend.text=element_text(size=11),
          legend.text.align=0,
          legend.title=element_text(face="bold")) +
        labs(x=paste(x_label, "for", query_title, "(Query)"), 
             y=paste(y_labels_query[i], " for", query_title, "(Query)"))
    }
    
    if(!is.null(match_min) & is.numeric(match_min)) {
      GG_Plot_List[[plotname]] <- GG_Plot_List[[plotname]] +
        geom_segment(data=DTW_DistRank[[i]][seq(1, nrow(DTW_DistRank[[i]]), sample_freq),], aes_string(x=colnames(DTW_DistRank[[i]])[3], y=colnames(DTW_DistRank[[i]])[5], xend=colnames(DTW_DistRank[[i]])[4], yend=colnames(DTW_DistRank[[i]])[6]), colour="black", alpha=match_vis, linetype=1)
    }
    
    
  }
  
  #Converting x_label to a version compatible with "parse"
  x_3way <- gsub(" ", "~", x_label)
  
  if(export_plots=="png") {
    for(i in names(GG_Plot_List)) {
      png(filename=paste0(export_path,"/", i,"_", format(Sys.time(), "%Y-%m-%d %Hhr %Mmin %Ssec"), ".png"),
          type="cairo", units="in",
          width=width, height=height,
          pointsize=point_size,
          res=dpi)
      print(GG_Plot_List[[i]])
      dev.off()
    }
    
    for(i in seq_along(alignment_list)) {
      png(filename=paste0(export_path,"/", names(GG_Plot_List)[i],"_3WAY_", format(Sys.time(), "%Y-%m-%d %Hhr %Mmin %Ssec"), ".png"),
          type="cairo", units="in",
          width=width, height=height,
          pointsize=point_size,
          res=dpi)
      
      Threeway_Alt(d=alignment_list[[i]], colours=colours[1:4],
                   xts = que_data[[query_subset[i]]], yts = ref_data[[ref_subset[i]]],
                   time_x = que_data[[x_var]], time_y = ref_data[[x_var]],
                   xlab = parse(text=paste0(x_3way,"~\"of\"~", y_labels_query[i], "~\"for\"~\"core\"~\"", query_title, "\"~\"(Query)\"")), 
                   ylab = parse(text=paste0(x_3way,"~\"of\"~", y_labels_ref[i], "~\"for\"~\"core\"~\"", ref_title, "\"~\"(Reference)\"")),
                   match.indices = match_subset, top_match = as.numeric(DTW_DistRank[[i]][,10])) #, top_match = as.numeric(DTW_DistRank[[i]][,10])
      dev.off()
    }
    
  } else if(export_plots=="pdf") {
    pdf(paste0(export_path,"/DTW_Alignment_Plots","_", format(Sys.time(), "%Y-%m-%d %Hhr %Mmin %Ssec"),".pdf"),
        width=width, height=height, pointsize=point_size)
    print(GG_Plot_List)
    
    for(i in seq_along(alignment_list)) {
      Threeway_Alt(d=alignment_list[[i]], colours=colours[1:4],
                   xts = que_data[[query_subset[i]]], yts = ref_data[[ref_subset[i]]],
                   time_x = que_data[[x_var]], time_y = ref_data[[x_var]],
                   xlab = parse(text=paste0(x_3way,"~\"of\"~", y_labels_query[i], "~\"for\"~\"core\"~\"", query_title, "\"~\"(Query)\"")), 
                   ylab = parse(text=paste0(x_3way,"~\"of\"~", y_labels_ref[i], "~\"for\"~\"core\"~\"", ref_title, "\"~\"(Reference)\"")),
                   match.indices = match_subset, top_match = as.numeric(DTW_DistRank[[i]][,10])) #, top_match = as.numeric(DTW_DistRank[[i]][,10])
    }
    
    dev.off()
  }
  
  print(GG_Plot_List)
  
  for(i in seq_along(alignment_list)) {
    Threeway_Alt(d=alignment_list[[i]], colours=colours[1:4],
                 xts = que_data[[query_subset[i]]], yts = ref_data[[ref_subset[i]]],
                 time_x = que_data[[x_var]], time_y = ref_data[[x_var]],
                 xlab = parse(text=paste0(x_3way,"~\"of\"~", y_labels_query[i], "~\"for\"~\"core\"~\"", query_title, "\"~\"(Query)\"")), 
                 ylab = parse(text=paste0(x_3way,"~\"of\"~", y_labels_ref[i], "~\"for\"~\"core\"~\"", ref_title, "\"~\"(Reference)\"")),
                 match.indices = match_subset, top_match = as.numeric(DTW_DistRank[[i]][,10])) #, top_match = DTW_DistRank[[i]][,10]
  }
  
  #return(list(que_data, ref_data))
  return(list(DTW_Data_List, alignment_list,GG_Plot_List))
}
