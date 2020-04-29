library(reticulate)
library(SingleCellExperiment)
library(tensorflow)
library(cellassign)

# Set up reticulate connection
use_python("/Users/ricard/anaconda3/envs/gpflow_v2/bin/python", required = T)
py_config()

#####################
## Define settings ##
#####################

# Load default settings
source("/Users/ricard/gastrulation10x/settings.R")

# Define I/O
io$outdir <- paste0(io$basedir,"/results/cellassign")

# Define options

# Testing mode
opts$test <- TRUE

# Cell types to use
opts$groups <- c(
  "Epiblast",
  "Primitive Streak",
  "ExE ectoderm",
  "Visceral endoderm",
  "ExE endoderm",
  "Nascent mesoderm",
  "Neurectoderm",
  "Blood progenitors",
  "Mixed mesoderm",
  "ExE mesoderm",
  "Pharyngeal mesoderm",
  "Caudal epiblast",
  "PGC",
  "Mesenchyme",
  "Haematoendothelial progenitors",
  "Surface ectoderm",
  "Gut",
  "Paraxial mesoderm",
  "Notochord",
  "Somitic mesoderm",
  "Caudal Mesoderm",
  "Erythroid",
  "Def. endoderm",
  "Parietal endoderm",
  "Allantois",
  "Anterior Primitive Streak",
  "Endothelium",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal cord",
  "Cardiomyocytes",
  "NMP",
  "Neural crest"
)

# Update metadata
sample_metadata %>% 
  setnames("celltype2","group")

# Acitvate test mode
if (isTRUE(opts$test)) {
  print("Testing mode, subsetting cells...")
  opts$groups <- head(opts$groups,n=3)
  sample_metadata <- sample_metadata %>% 
    .[group%in%opts$groups] %>% split(.,.$group) %>% 
    map(~ head(.,n=100)) %>% rbindlist
}
table(sample_metadata$group)


###############
## Load data ##
###############

# Load RNA expression data
sce <- readRDS(io$rna.sce)[,sample_metadata$cell]
sce$group <- sample_metadata$group

# Load marker genes
marker_genes.dt <- fread(io$marker_genes.lenient) %>%
  .[group%in%opts$groups]


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

# Covariate matrix (TO-DO: ADD BATCH INFORMATION)

# Run
fit <- cellassign(sce, 
  marker_gene_info = bmat, 
  min_delta = 1,
  # X = X,
  s = s, 
  # learning_rate = 1e-2, 
  shrinkage = TRUE,
  verbose = FALSE
)
print(fit)

##################
## Query output ##
##################

# celltypes(fit, assign_prob = 0.95)

# Plot heatmap of cell type probabilities
pdf(sprintf("%s/heatmap_probabilities.pdf",io$outdir), width = 8, height = 8)
pheatmap::pheatmap(cellprobs(fit))
dev.off()

# Compare to ground truth
print(table(sce$group, celltypes(fit)))

##########
## Save ##
##########

saveRDS(fit, paste0(io$outdir,"/cellassign_fit.rds"))