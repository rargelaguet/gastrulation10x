###################
## Load metadata ##
###################

# opts$celltypes <- opts$celltypes %>% stringr::str_replace_all("/", "_") %>% stringr::str_replace_all("_", " ") 

# Update metadata
sample_metadata <- sample_metadata %>%
  .[,celltype:=stringr::str_replace_all(celltype,"/", "_")] %>%
  .[,celltype:=stringr::str_replace_all(celltype,"_", " ")] %>%
  # .[celltype%in%opts$celltypes] %>%
  .[stage%in%args$stages]

metadata_query <- sample_metadata %>% .[sample%in%args$test_samples]
metadata_atlas <- sample_metadata %>% .[!sample%in%args$test_samples]
  
# Subset cells
if (isTRUE(opts$test_mode)) {
  metadata_query <- metadata_query %>% split(.$celltype) %>% map(~ head(.,n=50)) %>% rbindlist
  metadata_atlas <- metadata_atlas %>% split(.$celltype) %>% map(~ tail(.,n=50)) %>% rbindlist
}

###############
## Load data ##
###############

sce <- readRDS(io$rna.sce)

# remove non-expressed genes
sce <- sce[rowMeans(logcounts(sce))>0,]

# Define query SingleCellExperiment
sce.query <- sce[,metadata_query$cell]
colData(sce.query) <- metadata_query %>% tibble::column_to_rownames("cell") %>% DataFrame

# Define atlas SingleCellExperiment
sce.atlas <- sce[,metadata_atlas$cell]
colData(sce.atlas) <- metadata_atlas %>% tibble::column_to_rownames("cell") %>% DataFrame

#########################
## Subset marker genes ##
#########################

# Load gene markers to be used as HVGs
# marker_genes.dt <- fread(io$marker_genes)
# marker_genes.dt <- marker_genes.dt[,head(.SD,n=50),by="celltype"]
# marker_genes <- unique(marker_genes.dt$ens_id)
# marker_genes <- marker_genes[marker_genes%in%genes.intersect]
# 
# stopifnot(all(marker_genes%in%rownames(sce_atlas)))
# stopifnot(all(marker_genes%in%rownames(sce_query)))
# 
# # Update SingleCellExperiment objects
# sce_query <- sce_query[marker_genes,]
# sce_atlas <- sce_atlas[marker_genes,]
