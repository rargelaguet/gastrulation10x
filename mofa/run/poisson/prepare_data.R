library(Seurat)
library(data.table)
library(purrr)

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
io$outdir <- paste0(io$basedir,"/data/counts")


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

# opts$embryos <- c("E6.5_1", "E6.5_5", "E7.0_10", "E7.0_14", "E7.25_26", "E7.25_27")
opts$embryos <- c("E7.0_10")

opts$hvg <- 2500

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

# Subset lineages
cells.use <- meta.data %>%   
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>%
  .[stage_lineage%in%opts$stage_lineage,cell]
seurat <- subset(seurat, cells=cells.use)

# Subset stages
seurat@meta.data$stage_sample <- paste(seurat@meta.data$stage,seurat@meta.data$sample, sep="_")
cells <- rownames(seurat@meta.data[seurat@meta.data$stage_sample %in% opts$embryos,])
seurat <- subset(seurat, cells=cells)

###############
## Normalize ##
###############

seurat <- NormalizeData(seurat)

##################################
## Select highly variable genes ##
##################################

seurat <- FindVariableFeatures(seurat, nfeatures=opts$hvg, 
  selection.method = "mvp",
  mean.cutoff = c(0.1,2.5)
)
VariableFeaturePlot(seurat)

# seurat <- subset(seurat, features=seurat@assays$RNA@var.features)
hvg <- seurat@assays$RNA@var.features

###################
## Create matrix ##
###################

seurat@meta.data$stage_sample <- paste(seurat@meta.data$stage,seurat@meta.data$sample, sep="_")

m_counts <- as.matrix( seurat@assays$RNA@counts[hvg,] )
m_logcounts <- as.matrix( seurat@assays$RNA@data[hvg,] )

# Split by embryo
# m_counts <- split(as.data.frame(Matrix::t(seurat@assays$RNA@counts[hvg,])), seurat@meta.data$stage_sample)
# m_logcounts <- split(as.data.frame(Matrix::t(seurat@assays$RNA@data[hvg,])), seurat@meta.data$stage_sample)

# for (i in 1:length(m_counts)) {
#   m_counts[[i]]$sample <- NULL; m_counts[[i]] <- as.matrix(m_counts[[i]])
#   # m_logcounts[[i]]$sample <- NULL; m_logcounts[[i]] <- as.matrix(m_logcounts[[i]])
# }

##########
## Save ##
##########

# for (i in names(m_counts)) {
# 
#   # counts
#   filename <- sprintf("%s/%s_counts.txt",io$outdir,i)
#   write.table(m_counts[[i]], filename, row.names=T, col.names=T, quote=F, sep="\t")
# 
#   # logcounts
#   filename <- sprintf("%s/%s_logcounts.txt",io$outdir,i)
#   write.table(m_logcounts[[i]], filename, row.names=T, col.names=T, quote=F, sep="\t")
# }
# 
# system(sprintf("gzip -f %s/*.txt",io$outdir))
