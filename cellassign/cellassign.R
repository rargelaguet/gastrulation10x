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
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(cellassign))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--stages',    type="character",   nargs='+',  help='Stage(s)')
p$add_argument('--outdir',    type="character",               help='Output directory')
p$add_argument('--test',      action = "store_true",          help='Testing mode')
args <- p$parse_args(commandArgs(TRUE))

## START TEST 
# args$stages <- c("E6.5", "E6.75")
# args$outdir <- "/Users/ricard/data/gastrulation10x/results/cellassign"
# args$test <- TRUE
## END TEST

#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation10x/settings.R")
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/gastrulation10x/settings.R")
}  

# I/O
dir.create(args$outdir, show); dir.create(paste0(args$outdir,"/pdf"))

# Cell types to use
opts$celltypes <- opts$celltypes.1

# Maximum of M marker genes per cell type (sorted according to marker score)
opts$max.genes <- 50

# Update metadata
sample_metadata <- sample_metadata %>% 
  .[stage%in%args$stages] %>%
  setnames("celltype","group")


# Acitvate test mode
if (isTRUE(args$test)) {
  print("Testing mode, subsetting cells...")
  opts$celltypes <- sample(opts$celltypes,10)
  sample_metadata <- sample_metadata %>% 
    .[group%in%opts$celltypes] %>% split(.,.$group) %>% 
    map(~ head(.,n=50)) %>% rbindlist
}

# Minimum number of cells per lineage
opts$min.cells <- 50
foo <- table(sample_metadata$group)
if (any(foo<opts$min.cells)) {
  warning(sprintf("Some lineages have less than %d cells, removing them...",opts$min.cells))
  opts$celltypes <- names(which(foo>=opts$min.cells))
}
sample_metadata <- sample_metadata[group%in%opts$celltypes]


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
print(nrow(sce))

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
pdf(sprintf("%s/pdf/heatmap_probabilities.pdf",args$outdir), width = 8, height = 8)
pheatmap::pheatmap(cellprobs(fit))
dev.off()

# Compare to ground truth
foo <- table(sce$group, celltypes(fit))
pdf(sprintf("%s/pdf/heatmap_assignments.pdf",args$outdir), width = 8, height = 8)
pheatmap::pheatmap(foo, cluster_rows = F, cluster_cols = F)
dev.off()

##########
## Save ##
##########

outfile <- sprintf("%s/cellassign_fit_%s.rds",args$outdir,paste(args$stages,collapse="-"))
print(outfile)
saveRDS(fit, outfile)