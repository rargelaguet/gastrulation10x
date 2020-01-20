library(Seurat)
library(scran)
library(SingleCellExperiment)
library(data.table)
library(purrr)
library(umap)

matrix.please <- function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}

################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation10x/mofa/load_settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation10x/mofa/load_settings.R")
} else {
  stop("Computer not recognised")
}

io$seurat <- sprintf("%s/processed/seurat_E6.5_to_E7.5.rds",io$basedir)
io$outdir <- paste0(io$basedir,"/mofa")
# io$outfile <- paste0(io$outdir,"/data/E6.5-E7.0.txt")


####################
## Define options ##
####################

# Define stage and lineages
opts$stage_lineage10x_2 <- c(
	"E6.5_Epiblast",
  "E6.5_ExE ectoderm",
  "E6.5_ExE endoderm",
  "E6.5_Primitive Streak",
  "E6.5_Mesoderm",

  "E6.75_Epiblast",
  "E6.75_ExE ectoderm",
  "E6.75_ExE endoderm",
  "E6.75_Primitive Streak",
  "E6.75_Mesoderm",

  "E7.0_Epiblast",
	"E7.0_Mesoderm",
	"E7.0_ExE endoderm",
  "E7.0_ExE ectoderm",
  "E7.0_Primitive Streak",
	
	"E7.25_Epiblast",
	"E7.25_Mesoderm",
	"E7.25_ExE endoderm",
	"E7.25_ExE ectoderm",
	"E7.25_Primitive Streak"
)

###############
## Load data ##
###############

# Load seurat object
seurat <- readRDS(io$seurat)

# Load cell metadata
meta.data <- fread(io$sample.metadata)

#################
## Filter data ##
#################

# Filter out small lineages
# small_lineages <- names(which(table(seurat@meta.data$celltype) < 5))
# cells.use <- rownames(seurat@meta.data[!seurat@meta.data$celltype %in% small_lineages,])

# Subset preselected lineages
cells.use <- meta.data %>%   
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>%
  .[stage_lineage%in%opts$stage_lineage,cell]
seurat <- subset(seurat, cells=cells.use)

# Subset stages
seurat@meta.data$stage_sample <- paste(seurat@meta.data$stage,seurat@meta.data$sample, sep="_")
cells <- rownames(seurat@meta.data[seurat@meta.data$stage_sample %in% c("E6.5_1", "E6.5_5", "E7.0_10", "E7.0_14", "E7.25_26", "E7.25_27"),])
seurat <- subset(seurat, cells=cells)

#########################
## Normalise and scale ##
#########################

# Normalise
seurat <- NormalizeData(seurat, scale.factor=1000)

# Scale and regress out technical effects (batches)
seurat <- ScaleData(
  seurat,
  vars.to.regress = c("sample"),
  do.scale = FALSE, do.center = FALSE
)

##############################
## Dimensionality reduction ##
##############################

# tmp <- FindVariableFeatures(seurat, nfeatures=5000)
# tmp <- subset(tmp, features=tmp@assays$RNA@var.features)
# tmp <- ScaleData(tmp, do.scale = FALSE, do.center = FALSE)
# tmp <- ScaleData(tmp, do.scale = TRUE, do.center = FALSE, vars.to.regress = "sample")
# tmp <- Seurat::RunPCA(tmp)
# tmp <- Seurat::RunUMAP(tmp, reduction="pca", dims=1:15)

# Idents(tmp) <- tmp@meta.data$sample
# Idents(tmp) <- tmp@meta.data$celltype
# Seurat::DimPlot(tmp, reduction="pca", dims = c(5,6))
# Seurat::DimPlot(tmp, reduction="umap", dims = c(1,2))


##################################
## Select highly variable genes ##
##################################

seurat <- FindVariableFeatures(seurat, nfeatures=5000)
seurat <- subset(seurat, features=seurat@assays$RNA@var.features)


###################
## Create matrix ##
###################

# Split by batch
# tmp <- as.data.frame( t(seurat@assays$RNA@scale.data) )
tmp <- as.data.frame( t(seurat@assays$RNA@data) )
seurat@meta.data$stage_sample <- paste(seurat@meta.data$stage,seurat@meta.data$sample, sep="_")
m_batch <- split(tmp, seurat@meta.data$stage_sample)
for (i in 1:length(m_batch)) {
  m_batch[[i]]$sample <- NULL
  m_batch[[i]] <- as.matrix(m_batch[[i]])
}

lapply(m_batch,dim)


##########
## Save ##
##########

# Save seurat object
# save(seurat, file=io$outfile.seurat)

# Save matrices
for (i in names(m_batch)) {
  filename <- sprintf("%s/data/%s.txt",io$outdir,i)
  write.table(m_batch[[i]], filename, row.names=T, col.names=T, quote=F, sep="\t")
  system(sprintf("gzip -f %s",filename))
}

