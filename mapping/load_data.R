###################
## Load metadata ##
###################

# opts$celltypes <- opts$celltypes %>% stringr::str_replace_all("/", "_") %>% stringr::str_replace_all("_", " ") 

# Update metadata
sample_metadata <- sample_metadata %>% .[stage%in%args$stages]
if (!is.null(args$celltypes)) {
  sample_metadata <- sample_metadata %>%.[celltype%in%args$celltypes]
}

metadata_query <- sample_metadata %>% .[sample%in%as.character(args$test_samples)]
metadata_atlas <- sample_metadata %>% .[!sample%in%as.character(args$test_samples)]
  
# Subset cells
if (isTRUE(args$test)) {
  metadata_query <- metadata_query %>% split(.$celltype) %>% map(~ head(.,n=50)) %>% rbindlist
  metadata_atlas <- metadata_atlas %>% split(.$celltype) %>% map(~ tail(.,n=50)) %>% rbindlist
}

###############
## Load data ##
###############

# Load SingleCellExperiment
sce <- readRDS(io$rna.sce)

# remove lowly expressed genes
sce <- sce[rowSums(counts(sce))>10,]

# Define query
sce.query <- sce[,metadata_query$cell]
colData(sce.query) <- metadata_query %>% tibble::column_to_rownames("cell") %>% DataFrame

# Define atlas
sce.atlas <- sce[,metadata_atlas$cell]
colData(sce.atlas) <- metadata_atlas %>% tibble::column_to_rownames("cell") %>% DataFrame

rm(sce)

#################
## Subset HVGs ##
#################

# Select HVGs
# hvg <- getHVGs(sce.atlas, block=as.factor(sce.atlas$sample), p.value = 0.10)
# sce.query <- sce.query[hvg,]
# sce.atlas <- sce.atlas[hvg,]

#########################
## Subset marker genes ##
#########################

# # Load gene markers to be used as HVGs
# marker_genes.dt <- fread(io$marker_genes)
# marker_genes.dt <- marker_genes.dt[,head(.SD,n=50),by="celltype"]
# marker_genes <- unique(marker_genes.dt$ens_id)
# 
# stopifnot(all(marker_genes%in%rownames(sce.atlas)))
# stopifnot(all(marker_genes%in%rownames(sce.query)))
# 
# # Update SingleCellExperiment objects
# sce.query <- sce.query[marker_genes,]
# sce.atlas <- sce.atlas[marker_genes,]
