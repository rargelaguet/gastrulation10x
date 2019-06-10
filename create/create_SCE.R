library(data.table)
library(purrr)
library(scater)
library(Seurat)
# library(DropletUtils)
# library(scran)

read10xCounts.rna <- function (directory) {
  
  barcode.loc <- file.path(directory, "meta.tab")
  gene.loc <- file.path(directory, "genes.tsv")
  matrix.loc <- file.path(directory, "raw_counts.mtx")
  
  data_mat <- Matrix::readMM(matrix.loc) %>% as.matrix
  cell.info <- read.delim(barcode.loc, header = TRUE, colClasses = "character", stringsAsFactors = FALSE, sep="\t") %>% DataFrame()
  gene.info <- read.delim(gene.loc, header = FALSE, colClasses = "character", stringsAsFactors = FALSE, sep="\t") %>% DataFrame()
  colnames(gene.info) <- c("ens_id","gene_id")
  
  rownames(data_mat) <- gene.info$gene_id
  colnames(data_mat) <- cell.info$sample
  SingleCellExperiment(list(counts = data_mat), rowData = gene.info, colData = cell.info)
}

# Load data
sce <- read10xCounts.rna("/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/atlas")
dim(sce)

# Load size factors
# sf <- read.table("/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/atlas/sizefactors.tab")
# sizeFactors(sce) <- sf; sce$sizeFactor <- sf
  
############
## Filter ##
############

# no need to remove anything EXCEPT for the meta$doublet and meta$stripped indicated cells

## RNA ##
# Filter cell type:
sce <- sce[,sce$doublet==F & sce$stripped==F]

saveRDS(sce, "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/atlas/sce.rds")

###################
## Create Seurat ##
###################

sce <- readRDS("/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/atlas/sce.rds")
srat <- as.seurat(sce)

srat <- NormalizeData(srat, scale.factor=1000)

srat <- ScaleData(
  srat, 
  vars.to.regress="nUMI",
  model.use = "linear",
  do.scale = FALSE,
  do.center = TRUE, 
  do.par = TRUE, num.cores = 4
)
saveRDS(srat, "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/atlas/seurat.rds")

stop()

#############
## Process ##
#############
sce <- readRDS("/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/atlas/sce.rds")

sce.norm <- normalize(sce, exprs_values="counts")

sce = calculateQCMetrics(sce)

# size_factors = computeSumFactors(sce, positive=TRUE, sf.out=T)



saveRDS("/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/atlas/sce.rds")

##############################
## Dimensionality reduction ##
##############################

plotPCA(sce, colour_by="group", exprs_values="counts", ncomponents=5)

# plotTSNE(sce_rna.car_filt, colour_by="group", exprs_values="counts")

# foo = estimateSizeFactors(sce)
# cds = estimateDispersions(cds)
# cds = detectGenes(cds, 0.1)

# cds = reduceDimension(
#   cds, max_components = 2, num_dim = 15, norm_method = "log",
#   reduction_method = 'tSNE', 
#   residualModelFormulaStr = "~ experiment + num_genes_expressed", verbose = T)


plotPCA(sce, colour_by="group", ncomponents=5, run_args=list(exprs_values="logcounts"))
plotTSNE(sce, colour_by="group", add_ticks=F, run_args=list(exprs_values="logcounts", initial_dims=10))
# plotUMAP(sce, colour_by="group", add_ticks=F, run_args=list(exprs_values="counts"))

### TRY SEURAT ###

library(Seurat)

srat <- as.seurat(sce)

# srat <- CreateSeuratObject(
#   raw.data = matrix_filt, 
#   min.genes = 0, 
#   min.cells = 0,
#   project = "test", 
#   meta.data = metadata_filt %>% as.data.frame %>% tibble::column_to_rownames("sample"), 
#   normalization.method = "LogNormalize", 
#   do.scale=FALSE, 
#   do.center=FALSE
# )
srat <- NormalizeData(srat, scale.factor=1000)

srat@data[1:5,1:5]

srat <- ScaleData(
  srat, 
  vars.to.regress="total_counts",
  model.use = "linear",
  do.scale = FALSE,
  do.center = TRUE, 
  scale.max = 10
)

srat@raw.data[10:15,10:15]
srat@scale.data[1:5,1:5]
srat@data <- srat@scale.data
hist(colSums(srat@data))

foo <- tail(sort(apply(srat@scale.data,1,var)), n=1000)

Seurat::MeanVarPlot(srat)
srat <- Seurat::RunPCA(srat, pc.genes=names(foo), pcs.compute = 15)
srat <- Seurat::RunTSNE(srat)
srat <- Seurat::RunUMAP(srat)
Seurat::DimPlot(srat, reduction.use="pca", group.by = "group", dim.1 = 1, dim.2 = 2)
Seurat::DimPlot(srat, reduction.use="tsne", group.by = "group")
TSNEPlot(srat, group.by="group")
# Seurat::RunUMAP()

############
##########
## Save ##
##########

saveRDS(sce.rna.combined_filt, "/Users/ricard/data/Cao2018/HEK293T_NIH373_A549/combined/rna/parsed/sce.rds")
saveRDS(sce.atac.combined_filt, "/Users/ricard/data/Cao2018/HEK293T_NIH373_A549/combined/atac/parsed/sce.rds")