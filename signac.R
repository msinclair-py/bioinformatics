## create_object
## creates a new signac object
##
## arguments:
## * counts = matrix of raw counts with genes in the rows and cells in the columns
## * min.cells = remove genes that were detected in less than min.cells cells
## * min.features = remove cells that contain less than min.features detected genes
##
## returns values: signac object, implemented as a list composed of:
## * assay_data = counts matrix after filtering genes and cells
## * meta_data = a dataframe storing metadata for each of the cells
## * reductions = a list of dataframe that store coordinates from different dimension reduction methods

create_object = function(counts, min.cells = 0, min.features = 0) {
  
  ## genes to keep
  keep_genes = rowSums(counts > 0) >= min.cells
  
  ## cells to keep
  keep_cells = colSums(counts > 0) >= min.features
  
  ## create components of object
  assay_data = counts[keep_genes, keep_cells]
  meta_data = data.frame(nCount_RNA = colSums(assay_data))
  reductions = list(PCA = data.frame(X1 = rep(NA, nrow(meta_data))),
                    UMAP = data.frame(X1 = rep(NA, nrow(meta_data))))
  
  ## return object
  obj = list(assay_data = assay_data,
             meta_data =  meta_data,
             reductions = reductions)
  return(obj)
}

## percentage_feature_set
## creates percentage of counts from a set of features
##
## arguments:
## * obj = a signac object
## * features = a vector of row numbers of the features belonging to the set of interest
## 
## return values:
## * Percentage of counts from each cell coming from the features in features

percentage_feature_set = function(obj, features) {
  percent = colSums(obj$assay_data[features,]) / colSums(obj$assay_data)
  return(percent)
}

## subset_obj
## subsets cells from a signac object
##
## arguments:
## * obj = a signac object
## * subset = a vector of column numbers belonging to the cells that should be kept
##
## return values: signac object with cells not in subset removed from signac object's
## * assay_data
## * reductions$PCA and reductions$UMAP
## * meta_data

subset_obj = function(obj, subset) {
  
  ## remove cells
  obj$assay_data = obj$assay_data[, subset]
  obj$meta_data = obj$meta_data[subset,]
  obj$reductions$PCA = obj$reduction$PCA[subset,, drop = FALSE]
  obj$reductions$UMAP = obj$reductions$UMAP[subset,, drop = FALSE]
  
  ## return object
  return(obj)
}

## normalize_data
## normalizes expression data across cells by:
## 1. converting to transcripts per scale.factor
## 2. natural log-transforming after the addition of 1
##
## arguments:
## * obj = a signac object
## * scale.factor = a scaling factor
## 
## return values: signac object with counts in the assay_data matrix

normalize_data = function(obj, scale.factor = 10000) {
  
  ## divide each row of assay_data by the total number of RNA counts from that cell
  ## scale by scale.factor
  ## taking the transpose allows us to use vector division, which is fast in R
  obj$assay_data = t(t(obj$assay_data) / obj$meta_data$nCount_RNA * scale.factor)
  
  ## add 1 and log-transform
  obj$assay_data = log(obj$assay_data + 1)
  
  ## return object
  return(obj)
}

## regress_out
## replaces expression data with residuals after regressing out unwanted variation
##
## arguments:
## * obj = a signac object
## * vars.to.regress = a vector of variable names stored in the meta_data that represent unwanted sources of variation
##
## return values: signac object with assay_data replaced by residuals

regress_out = function(obj, vars.to.regress) {
  
  meta_data_vars = data.frame(obj$meta_data[, vars.to.regress])
  
  ## for each gene:
  for (gene in 1:nrow(obj$assay_data)) {
    
    ## run regression
    fit = lm(obj$assay_data[gene,] ~ ., data = meta_data_vars)
    ## calculate residuals and replace assay_data
    obj$assay_data[gene,] = residuals(fit)
  }
  
  ## return object
  return(obj)
}

## find_variable_features
## find genes with highest IQR
##
## arguments:
## * obj = a signac object
## * nfeatures = number of most variable features to find; default value = 2000
##
## return values: signac object with new element variable_genes containing the names of most variable genes

find_variable_features = function(obj, nfeatures = 2000) {
  
  ## calculate IQRs
  IQRs = apply(obj$assay_data, 1, IQR)
  ## find nfeatures genes with largest IQRs
  genes = rownames(obj$assay_data)[order(-IQRs)[1:nfeatures]]
  ## update meta_data
  obj$variable_genes = genes
  
  return(obj)
}

## run_pca
## runs PCA on centered and scaled assay_data
##
## arguments:
## * obj = a signac object
## * npcs = number of principal components to keep
##
## return values: signac object with updated reductions$PCA

run_pca = function(obj, npcs) {
  
  ## calculate PCs
  obj$reductions$PCA = prcomp(t(obj$assay_data[obj$variable_genes,]),
                              center = TRUE,
                              scale. = TRUE)$x[, 1:npcs]
  
  return(obj)
}

## run_umap
## runs UMAP on reductions$PCA
##
## arguments:
## * obj = a signac object
##
## return values: signac object with updated reductions$UMAP

library(umap)

run_umap = function(obj) {
  
  ## calculate UMAP
  obj$reductions$UMAP = umap(obj$reductions$PCA)$layout
  return(obj)
  
}

## find_clusters
## finds gene clusters by:
## 1. constructing a shared nearest neighbor network using reductions$PCA
## 2. applying the louvain method
##
## arguments:
## * obj = a signac object
## * k = number of neighbors to consider to calculate the shared nearest neighbors
## * kt = two points are connected in the sNN if they share kt or more neighbors
##
## return values: signac object with cluster memberships added to meta_data

library(dbscan)
library(igraph)

find_clusters = function(obj, k, kt) {
  
  ## make snn
  snn = sNN(obj$reductions$PCA, k = k, kt = kt)
  ## make undirected graph from adjacency list
  g = graph_from_adj_list(adjacencylist(snn),
                          mode = "all",
                          duplicate = FALSE)
  ## update meta_data object
  obj$meta_data[["cluster"]] = as.factor(membership(cluster_louvain(g)))
  
  return(obj)
}

## feature_plot
## visualizes the expression levels of specific genes of interest
##
## arguments:
## * obj = a signac object
## * features = a vector of gene names to visualize
##
## return values: list of ggplot objects that visualize the expression of each gene on the umap plot

library(ggplot2)

feature_plot = function(obj, features) {
  
  features = intersect(features, rownames(obj$assay_data))
  plots = vector(mode = "list", length = length(features))
  
  for (i in 1:length(features)) {
    df = data.frame(UMAP1 = obj$reductions$UMAP[,1],
                    UMAP2 = obj$reductions$UMAP[,2],
                    Expression = obj$assay_data[features[i],])
    plots[[i]] = ggplot(data = df) +
      geom_point(mapping = aes(x = UMAP1, y = UMAP2, color = Expression))+
      ggtitle(features[i])
  }
  
  return(plots)
}

## find_markers
## identify genes that are significantly differentially expressed in a cluster of interest than in any other clusters
##
## arguments:
## * obj = a signac object
## * ident.1 = the cluster of interest
##
## return values: dataframe where each row is a gene and each column with the following columns
## * p_val = p-value
## * ave_expr_diff = average expression in the cluster minus average expression in all other cells
## * p_val_adj = Bonferroni-adjusted p-value
## the rows are ordered in ascending order by p-value

find_markers = function(obj, ident.1) {
  
  clust = obj$meta_data[["cluster"]] == ident.1
  
  p_val = apply(obj$assay_data, 1, function(x) {
    wilcox.test(x[clust], x[!clust])$p.value
  })
  p_val_adj = p.adjust(p_val, method = "bonferroni")
  
  ave_expr_diff =
    rowMeans(obj$assay_data[, clust]) -
    rowMeans(obj$assay_data[, !clust])
  
  results = data.frame(p_val = p_val,
                       ave_expr_diff = ave_expr_diff, 
                       p_val_adj = p_val_adj)
  
  return(results[order(results$p_val),])
}
