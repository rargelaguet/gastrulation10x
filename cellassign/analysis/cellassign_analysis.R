#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation10x/settings.R")
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/gastrulation10x/settings.R")
}  


# Define stages
opts$stages <- c("E6.5")

# Update metadata
sample_metadata <- sample_metadata %>% 
  .[stage%in%opts$stages] %>%
  setnames("celltype","group")

# (TO-DOF IX THIS) Minimum number of cells per lineage
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

# Load precomputed cellassign model
io$model <- sprintf("%s/cellassign_fit_%s.rds",io$cellassign.dir,paste(opts$stages,collapse="-"))
fit <- readRDS(io$model)
print(fit)

##################
## Query output ##
##################

# Plot heatmap of cell type probabilities
pdf(sprintf("%s/pdf/heatmap_probabilities_%s.pdf",args$outdir,paste(args$stages,collapse="-")), width = 8, height = 8)
pheatmap::pheatmap(cellprobs(fit))
dev.off()

# Compare to ground truth
foo <- prop.table(table(sample_metadata$group, celltypes(fit)), margin=1)
pdf(sprintf("%s/pdf/heatmap_assignments_%s.pdf",args$outdir,paste(args$stages,collapse="-")), width = 8, height = 8)
pheatmap::pheatmap(foo, cluster_rows = F, cluster_cols = F, main = "Rows = TRUE, Cols = PRED")
dev.off()

##########
## Save ##
##########


saveRDS(fit, outfile)