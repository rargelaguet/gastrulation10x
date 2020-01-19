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
  "E6.5_Mesoderm"

 #  "E6.75_Epiblast",
 #  "E6.75_ExE ectoderm",
 #  "E6.75_ExE endoderm",
 #  "E6.75_Primitive Streak",
 #  "E6.75_Mesoderm",

 #  "E7.0_Epiblast",
	# "E7.0_Mesoderm",
	# "E7.0_ExE endoderm",
 #  "E7.0_ExE ectoderm",
 #  "E7.0_Primitive Streak",
	
	# "E7.25_Epiblast",
	# "E7.25_Mesoderm",
	# "E7.25_ExE endoderm",
	# "E7.25_ExE ectoderm",
	# "E7.25_Primitive Streak"
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


##################
## Prepare data ##
##################


data <- rna_dt %>%
  # .[,stage_assay:=paste(stage,assay,sep="_")] %>%
  .[,c("id_rna","assay","gene","expr")] %>%
  setnames(c("sample","sample_group","feature","value")) %>%
  .[,c("feature_group"):="RNA"]

# fwrite(data,paste0(io$outdir,"/data.txt"), col.names=T, quote=F, sep="\t")
# system(sprintf("pigz -f %s",paste0(io$outdir,"/data.txt")))

# data[,.N,by=c("sample_group","feature","feature_group")] %>% View


######################
## Train MOFA model ##
######################

# Create input data.frame

# lapply(object@input_data, function(x) dim(x[[2]]) )

data_opts <- get_default_data_options(object)

model_opts <- get_default_model_options(object)
model_opts$num_factors <- 10
model_opts$spikeslab_z <- FALSE
model_opts$spikeslab_w <- TRUE
model_opts$ard_w <- TRUE
model_opts$ard_z <- FALSE


train_opts <- get_default_training_options(object)
train_opts$maxiter <- 10
train_opts$convergence_mode <- "medium"
train_opts$seed <- 1
train_opts$startELBO <- 100
train_opts$gpu_mode <- TRUE

object <- prepare_biofam(object, 
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)


model <- run_biofam(object, outfile=outfile)
