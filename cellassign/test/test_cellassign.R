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

# opts$stages <- c("E6.5")
opts$stages <- c(
  "E6.5",
  "E6.75",
  "E7.0",
  "E7.25",
  "E7.5",
  "E7.75",
  "E8.0",
  "E8.25",
  "E8.5",
  "mixed_gastrulation"
)

opts$test <- TRUE

# Cell types to use
# opts$celltypes <- opts$celltypes.1
opts$celltypes <- c(
  "Erythroid1", "Erythroid2",
  "Visceral endoderm", "ExE endoderm",
  "Epiblast", "Primitive Streak"
)

# Maximum of M marker genes per cell type (sorted according to marker score)
opts$max.genes <- 25

# Update metadata
sample_metadata <- sample_metadata %>% 
  .[,celltype:=stringr::str_replace_all(celltype,"_", " ")] %>%
  .[stage%in%opts$stages & celltype%in%opts$celltypes] %>%
  setnames("celltype","group")

# Subset cells
if (isTRUE(opts$test)) {
  sample_metadata <- sample_metadata %>% 
    split(.$group) %>% map(~ head(.,n=100)) %>% rbindlist
}

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
  .[,celltype:=stringr::str_replace_all(celltype,"_", " ")] %>%
  setnames("celltype","group") %>%
  .[group%in%opts$celltypes] %>%
  setorder(group,-score)

# Sanity check
stopifnot(all(unique(marker_genes.dt$group)%in%opts$celltypes))

# Maximum number of marker genes per cell type
marker_genes.dt <- marker_genes.dt %>% 
  .[,head(.SD,n=opts$max.genes),by="group"]

table(marker_genes.dt$group)

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
bmat <- bmat[,colnames(bmat)!="other"]

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
celltype.pred <- celltypes(fit, assign_prob = 0.50)
foo <- table(sce$group, celltype.pred)
pheatmap::pheatmap(foo, cluster_rows = F, cluster_cols = F)


foobar <- merge(sample_metadata[,c("cell","group")], data.table(cell=colnames(sce), celltype.pred=celltype.pred), by="cell")
mean(foobar$group==foobar$celltype.pred)
