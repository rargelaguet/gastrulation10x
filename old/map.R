
# Script to map the scNMT-seq data to the 10x mouse reference
"""
To-do:
- Show where our three E7.5 lineages match

# You input two matrices: hvg_map and hvg_reference. These should have the same genes (rows) in the same order. 
I use normalised counts - but for your data, perhaps you want to log transform and z-score - perhaps try. 
# You also need to input metaA, which is the meta data for the atlas/reference. 
"""

library(purrr)
library(data.table)
library(SingleCellExperiment)
source("/Users/ricard/data/gastrulation10x/scripts/core_functions.R")



# Load reference
# sce <- load_data(normalise = TRUE, hvg = TRUE, corrected = FALSE, remove_doublets = FALSE)
# cData <- meta %>% tibble::remove_rownames() %>% tibble::column_to_rownames("cell") %>% DataFrame
# colData(sce) <- cData
sce <- readRDS("/Users/ricard/data/gastrulation10x/data/parsed/sce_parsed.rds")

# Load scNMT-seq data
scnmt <- readRDS("/Users/ricard/data/gastrulation/rna/parsed/sceset_scNMT.rds")
scnmt$stage_lineage <- paste(scnmt$stage,scnmt$lineage,sep="_")

###########################################
## Subset stage and lineages of interest ##
###########################################

sce_filt <- list()
scnmt_filt <- list()

## Subset cols ##

# stage_lineage
scnmt_filt[["E7.5"]] <- scnmt[,scnmt$stage_lineage%in%c("E7.5_Mesoderm")]
sce_filt[["E7.5"]] <- sce[,colData(sce)$stage%in%c("E6.5","E7.5","E8.0")]


## Subset rows ##

# overlapping genes
genes <- intersect(rownames(sce),rownames(scnmt_filt[["E7.5"]]))
scnmt_filt[["E7.5"]] <- scnmt_filt[["E7.5"]][genes,]
sce_filt[["E7.5"]] <- sce_filt[["E7.5"]][genes,]

# highly variable genes
sce_filt[["E7.5"]] <- sce_filt[["E7.5"]][rownames(sce_filt[["E7.5"]]) %in% hvg.list$E7.5,]

rm(sce)
gc(reset=T)


#############
## Mapping ##
#############

# Extract expression matrix
hvg_map <- exprs(scnmt_filt[["E7.5"]])
hvg_reference <- exprs(sce_filt[["E7.5"]])

NN <- 10 # number of nearest neighbours

dist_correl <- sqrt(0.5*(1-cor(as.matrix(exprs(scnmt_filt[["E7.5"]])), as.matrix(exprs(sce_filt[["E7.5"]])), method="pearson")))
# dist_correl <- sqrt(0.5*((1-cor(as.matrix(hvg_map),as.matrix(hvg_reference),method="pearson"))))
# rownames(dist_correl) <- colnames(hvg_map)
# colnames(dist_correl) <- colnames(hvg_reference)
knn <- apply(dist_correl,1,function(x) names(x)[order(x)[1:NN]])


### 
map_cells <- colnames(scnmt_filt[["E7.5"]])[scnmt_filt[["E7.5"]]$lineage=="Mesoderm"]
reference_cells <- as.character(knn[,map_cells])
sce_filt[["E7.5"]]$tag <- colnames(sce_filt[["E7.5"]]) %in% reference_cells
sce_pca <- runPCA(sce_filt[["E7.5"]], ncomponents=3, scale_features = F)
plotPCA(sce_pca, colour_by="stage", shape_by="tag")
