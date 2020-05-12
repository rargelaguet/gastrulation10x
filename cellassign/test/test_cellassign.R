# Set up reticulate connection
suppressPackageStartupMessages(library(reticulate))
if (grepl("ricard",Sys.info()['nodename'])) {
  use_python("/Users/ricard/anaconda3/envs/gpflow_v2/bin/python", required = TRUE)
} else if(grepl("ebi",Sys.info()['nodename'])){
  use_python("/nfs/research1/stegle/users/ricard/conda-envs/gpflow2/bin/python", required = TRUE)
}  
py_config()

# Load libraries
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tensorflow))
suppressPackageStartupMessages(library(cellassign))

#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation10x/settings.R")
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/gastrulation10x/settings.R")
}  
io$outdir <- "/Users/ricard/data/gastrulation10x/results/cellassign/test"
dir.create(io$outdir, showWarnings = F)

opts$stages <- c("E6.5")
opts$test <- TRUE

# Cell types to use
opts$celltypes <- opts$celltypes.1

# Maximum of M marker genes per cell type (sorted according to marker score)
opts$max.genes <- 50

# Update metadata
sample_metadata <- sample_metadata %>% 
  .[stage%in%opts$stages] %>%
  .[celltype%in%c("Epiblast","Primitive_Streak")] %>%
  setnames("celltype","group")

# Subset cels
sample_metadata <- sample_metadata %>% split(.$group) %>% 
    map(~ head(.,n=100)) %>% rbindlist

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
  .[group%in%opts$celltypes] %>%
  setorder(group,-score)

marker_genes.dt <- fread("/Users/ricard/data/gastrulation10x/results/differential/Epiblast_vs_Primitive_Streak.txt.gz") %>%
  .[sig==T] %>%
  .[,group:=ifelse(logFC>0,"Primitive_Streak","Epiblast")] %>%
  .[,c("group","ens_id")]


# Sanity check
stopifnot(all(unique(marker_genes.dt$group)%in%opts$celltypes))

# Maximum number of marker genes per cell type
marker_genes.dt <- marker_genes.dt[,head(.SD,n=opts$max.genes),by="group"]

######################
## Print statistics ##
######################

print("Number of marker genes per cell type")
print(marker_genes.dt[,.N,by="group"])

print("Total number of marker genes")
print(length(unique(marker_genes.dt$ens_id)))

print("Number of cells per cell type:")
table(sample_metadata$group)

print("Total number of cells")
print(ncol(sce))

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
pheatmap::pheatmap(cellprobs(fit))

# Compare to ground truth
foo <- table(sce$group, celltypes(fit))
pheatmap::pheatmap(foo, cluster_rows = F, cluster_cols = F)
