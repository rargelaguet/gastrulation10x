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

## Define I/O ##

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation10x/mofa/load_settings.R")
  io$basedir <- "/Users/ricard/data/gastrulation10x"
  io$gene.metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
} else {
  source("/homes/ricard/gastrulation10x/mofa/load_settings.R")
  io$basedir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x"
}
io$outdir <- paste0(io$basedir,"/mofa")
# io$outfile <- paste0(io$outdir,"/data/E6.5-E7.0.txt")


## Define options ##

# Define which stage and lineages to look at 
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
# io$seurat <- sprintf("%s/processed/seurat_%s.rds",io$basedir,opts$stage)
io$seurat <- sprintf("%s/processed/seurat/seurat_E6.5_to_E7.5.rds",io$basedir)
srat <- readRDS(io$seurat)

# Load sample metadata
# meta.data <- srat@meta.data
meta.data <- fread(io$sample.metadata)

#################
## Filter data ##
#################

## Filter out cells ##

# Filter out small lineages
if (is.null(opts$stage_lineage10x_2)) {
  small_lineages <- names(which(table(srat@meta.data$celltype) < 5))
  cells.use <- rownames(srat@meta.data[!srat@meta.data$celltype %in% small_lineages,])
# Filter out preselected lineages
} else {
  cells.use <- meta.data %>%   
    .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>%
    .[stage_lineage%in%opts$stage_lineage,cell]
}
srat <- subset(srat, cells=cells.use)

####################
## Normalize data ##
####################

# Normalize
srat <- NormalizeData(srat, scale.factor=1000)

# Scale and regress out batch effects
# srat <- ScaleData(
#   srat, 
#   # vars.to.regress=c("sample"),
#   # model.use = "linear",
#   do.scale = FALSE, do.center = FALSE
# )

# Select highly variable genes
srat <- FindVariableFeatures(srat, nfeatures=5000)
srat <- subset(srat, features=srat@assays$RNA@var.features)


###################
## Create matrix ##
###################

# Split by batch
# tmp <- as.data.frame( t(srat@assays$RNA@scale.data) )
tmp <- as.data.frame( t(srat@assays$RNA@data) )
srat@meta.data$stage_sample <- paste(srat@meta.data$stage,srat@meta.data$sample, sep="_")
m_batch <- split(tmp, srat@meta.data$stage_sample)
for (i in 1:length(m_batch)) {
  m_batch[[i]]$sample <- NULL
  m_batch[[i]] <- as.matrix(m_batch[[i]])
}

lapply(m_batch,dim)

###################################
## Dimensionality reduction test ##
###################################

# # Do fast PCA with prcomp_irba
# Z.irlba <- irlba::prcomp_irlba(x = t(m), n = 50)$x

# # t-SNE
# # tsne <- Rtsne(Z.irlba, check_duplicates=FALSE, pca=FALSE, theta=0.5, dims=2, initial_dims=model@dimensions$K)
# # Z.out <- tsne$Y

# # UMAP
# umap.defaults$n_neighbors <- 20
# umap.defaults$min_dist <- 0.7
# umap.out <- umap(Z.irlba, config = umap.defaults)
# Z.out <- umap.out$layout

# to.plot <- Z.out %>% as.data.table %>% .[,cell:=colnames(m)] %>%
#   merge(meta.data, by="cell")

# p <- ggplot(to.plot, aes(x=V1, y=V2, color=`celltype`)) +
#   geom_point(alpha=0.7, size=1.5) +
#   labs(x="", y="") +
#   scale_color_manual(values=opts$colors) +
#   theme_classic()
#   # theme(
#   #   legend.position = "none"
#   # )
# # p

# pdf(paste0(io$outdir,"/celltype.pdf"))
# print(p)
# dev.off()

# p <- ggplot(to.plot, aes(x=V1, y=V2, color=sample)) +
#   geom_point(alpha=0.7, size=1.5) +
#   labs(x="", y="") +
#   theme_classic()
# # p

# pdf(paste0(io$outdir,"/batch.pdf"))
# print(p)
# dev.off()


##########
## Save ##
##########

# write.table(t(m), io$outfile, row.names=T, col.names=T, quote=F, sep="\t")
# system(sprintf("gzip -f %s",io$outfile))

for (i in names(m_batch)) {
  filename <- sprintf("%s/data/%s.txt",io$outdir,i)
  write.table(m_batch[[i]], filename, row.names=T, col.names=T, quote=F, sep="\t")
  system(sprintf("gzip -f %s",filename))
}

