
####################################################################
## Plot dimensionality reduction of EB cells mapped to the atlas ##
####################################################################

library(SingleCellExperiment)
library(data.table)
library(purrr)
library(ggplot2)

source("/Users/ricard/NMT-seq_EB+ESC/rna/mapping/settings.R")

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
  .[!celltype%in%c("PGC","Parietal endoderm")] %>%
  .[,aggregated_celltype:=stringr::str_replace_all(celltype,aggregated_celltypes)]

# Load precomputed dimensionality reduction coordinates
umap <- meta_atlas %>%
  setnames(c("umapX","umapY"),c("V1","V2"))

################
## Parse data ##
################

umap <- umap[!celltype=="ExE ectoderm"]

###################################
## Plot dimensionality reduction ##
###################################

plot_df <- umap

p <- ggplot(data=plot_df, mapping=aes(x=V1, y=V2)) +
  ggrastr::geom_point_rast(aes(colour=aggregated_celltype), size=0.25) +
  scale_color_manual(values=aggregated_celltype_colours) +
  
  # ggrastr::geom_point_rast(aes(colour=celltype), size=0.25) +
  # scale_color_manual(values=celltype_colours) +
  
  labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

pdf(paste0(io$outdir,"/pdf/umap_v2_2.pdf"), width=5, height=5, useDingbats = F)
print(p)
dev.off()
