here::i_am("scanpy/convert_SingleCellExperiment_to_anndata.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(reticulate))

#####################
## Define settings ##
#####################

# I/O
io$outfile <- file.path(io$basedir,"processed/anndata.h5ad")

# Options
opts$python_path <- "/Users/argelagr/opt/anaconda3/envs/main/bin/python"

# Reticulate
reticulate::use_python(args$python_path, required = TRUE)
sc <- import("scanpy")

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[stripped==FALSE & doublet==FALSE]

###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(io$sce, cells=sample_metadata$cell, normalise = FALSE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

#############################################
## Convert SingleCellExperiment to AnnData ##
#############################################

adata_sce <- sc$AnnData(
    X   = t(counts(sce)),
    obs = as.data.frame(colData(sce)),
    var = data.frame(gene=rownames(sce), row.names=rownames(sce))
)
# adata_sce$obsm$update(umap = reducedDim(sce, "umap"))

adata_sce

##########################
## Parse anndata object ##
##########################

# Add dimensionality reduction
corrected_pca.mtx <- readRDS("/Users/argelagr/data/gastrulation10x/results/dimensionality_reduction/corrected_pcas.rds")$all %>% .[rownames(adata_sce$obs),]
adata_sce$obsm$pca_corrected_mnn <- corrected_pca.mtx

umap.mtx <- sample_metadata[,c("cell","umapX","umapY")] %>% matrix.please %>% .[rownames(adata_sce$obs),]
adata_sce$obsm$precomputed_umap <- umap.mtx

# Add stage colors
adata_sce$uns$update(stage_colors = opts$stage.colors[sort(unique(as.character(adata_sce$obs$stage)))])
adata_sce$uns["stage_colors"]

# Add cell type colors
adata_sce$uns$update(celltype_colors = opts$celltype.colors[sort(unique(as.character(adata_sce$obs$celltype)))])
adata_sce$uns["celltype_colors"]

# colPalette_celltypes = [opts["celltype_colors"][i.replace(" ","_").replace("/","_")] for i in sorted(np.unique(adata.obs['celltype']))]
# adata.uns['celltype'] = colPalette_celltypes
# colPalette_stages = [opts["stage_colors"][i.replace(" ","_").replace("/","_")] for i in sorted(np.unique(adata.obs['stage']))]
# adata.uns['stage_colors'] = colPalette_stages

##########
## Save ##
##########

adata_sce$write_h5ad(io$outfile)
