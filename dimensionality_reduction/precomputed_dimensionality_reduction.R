################################################
## Plot pre-computed dimensionality reduction ##
################################################

library(SingleCellExperiment)
library(data.table)
library(purrr)
library(ggplot2)

source("/Users/ricard/gastrulation10x/dimensionality_reduction/utils.R")

#####################
## I/O and options ##
#####################

io <- list()
io$path2atlas <- "/Users/ricard/data/gastrulation10x"
io$outdir <- "/Users/ricard/data/gastrulation10x/dimensionality_reduction"

####################
## Load 10x atlas ##
####################

# Load atlas metadata
meta_atlas <- fread(paste0(io$path2atlas, "/sample_metadata.txt")) %>%
  .[stripped==F & doublet==F] %>%
  .[,aggregated_celltype:=stringr::str_replace_all(celltype,aggregated_celltypes)]

# Load precomputed dimensionality reduction coordinates
umap <- meta_atlas %>%
  setnames(c("umapX","umapY"),c("V1","V2"))

################
## Parse data ##
################

# Remove specific lineages
plot_df <- umap[!celltype%in%c("ExE_ectoderm","PGC", "Parietal_endoderm")]

###################################
## Plot dimensionality reduction ##
###################################

p <- ggplot(data=plot_df, mapping=aes(x=V1, y=V2)) +
  ggrastr::geom_point_rast(aes(colour=aggregated_celltype), size=0.25) +
  # ggrastr::geom_point_rast(aes(colour=celltype), size=0.25) +
  scale_color_manual(values=colors) +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )

pdf(paste0(io$outdir,"/umap_forCarine.pdf"), width=8.5, height=5, useDingbats = F)
print(p)
dev.off()
