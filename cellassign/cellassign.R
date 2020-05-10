# Set up reticulate connection
library(reticulate)
if (grepl("ricard",Sys.info()['nodename'])) {
  use_python("/Users/ricard/anaconda3/envs/gpflow_v2/bin/python", required = TRUE)
} else if(grepl("ebi",Sys.info()['nodename'])){
  use_python("/nfs/research1/stegle/users/ricard/conda-envs/gpflow2/bin/python", required = TRUE)
}  
py_config()

library(SingleCellExperiment)
library(tensorflow)
library(cellassign)

#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation10x/settings.R")
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/gastrulation10x/settings.R")
}  
# io$marker_genes <- paste0(io$basedir,"/results/cellassign/marker_genes.txt.gz")

# Define I/O
io$outdir <- paste0(io$basedir,"/results/cellassign")
# Define options

# Cell types to use
opts$groups <- opts$celltypes.1

# Stages to use
opts$stages <- c(
  # "E6.5",
  # "E6.75",
  # "E7.0"
  # "E7.25",
  "E7.5"
  # "E7.75",
  # "E8.0",
  # "E8.25",
  # "E8.5",
  # "mixed_gastrulation"
)

# Maximum of M marker genes per cell type (sorted according to marker score)
opts$max.genes <- 50

# Update metadata
sample_metadata <- sample_metadata %>% 
  .[stage%in%opts$stages] %>%
  setnames("celltype","group")

# Testing mode
opts$test <- FALSE

# Acitvate test mode
if (isTRUE(opts$test)) {
  print("Testing mode, subsetting cells...")
  opts$groups <- sample(opts$groups,10)
  sample_metadata <- sample_metadata %>% 
    .[group%in%opts$groups] %>% split(.,.$group) %>% 
    map(~ head(.,n=50)) %>% rbindlist
}

# Minimum number of cells per lineage
opts$min.cells <- 50
foo <- table(sample_metadata$group)
if (any(foo<opts$min.cells)) {
  warning(sprintf("Some lineages have less than %d cells, removing them...",opts$min.cells))
  opts$groups <- names(which(foo>=opts$min.cells))
}
sample_metadata <- sample_metadata[group%in%opts$groups]
table(sample_metadata$group)

###############
## Load data ##
###############

# Load RNA expression data
sce <- readRDS(io$rna.sce)[,sample_metadata$cell]
sce$group <- sample_metadata$group

#######################
## Load marker genes ##
#######################

marker_genes.dt <- fread(io$marker_genes) %>%
  setnames("celltype","group") %>%
  .[group%in%opts$groups] %>%
  setorder(group,-score)

# Sanity check
stopifnot(all(unique(marker_genes.dt$group)%in%opts$groups))

# Maximum number of marker genes per cell type
marker_genes.dt <- marker_genes.dt[,head(.SD,n=opts$max.genes),by="group"]
print(marker_genes.dt[,.N,by="group"])

################
## cellassign ##
################

# Create binary membership matrix
marker_gene_list <- split(marker_genes.dt, by="group") %>% map(~ .$ens_id)
bmat <- marker_list_to_mat(marker_gene_list)

# Subset SingleCellExperiment
sce <- sce[rownames(bmat),]

# Extract size factors
s <- sizeFactors(sce)

# Run
fit <- cellassign(sce, 
  marker_gene_info = bmat, 
  min_delta = 1,
  s = s, 
  # learning_rate = 1e-2, 
  shrinkage = TRUE,
  verbose = FALSE
)
print(fit)

##################
## Query output ##
##################

# Plot heatmap of cell type probabilities
pdf(sprintf("%s/heatmap_probabilities.pdf",io$outdir), width = 8, height = 8)
pheatmap::pheatmap(cellprobs(fit))
dev.off()

# Compare to ground truth
foo <- table(sce$group, celltypes(fit))
pdf(sprintf("%s/heatmap_assignments.pdf",io$outdir), width = 8, height = 8)
pheatmap::pheatmap(foo, cluster_rows = F, cluster_cols = F)
dev.off()

##########
## Save ##
##########

outfile <- sprintf("%s/cellassign_fit_%s.rds",io$outdir,paste(opts$stages,collapse="-"))
print(outfile)
saveRDS(fit, outfile)