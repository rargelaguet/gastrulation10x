here::i_am("dimensionality_reduction/precomputed_dimensionality_reduction.R")

source(here::here("settings.R"))

#####################
## I/O and options ##
#####################

io$outdir <- paste0(io$basedir,"/results/dimensionality_reduction/pdf")

opts$aggregated.celltypes <- c(
  "Erythroid1" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid3" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors",
  "Rostral_neurectoderm" = "Neurectoderm",
  "Caudal_neurectoderm" = "Neurectoderm",
  "Anterior_Primitive_Streak" = "Primitive_Streak"
)

opts$stages <- c(
  "E6.5",
  "E6.75",
  "E7.0",
  "E7.25",
  "E7.5",
  "E7.75",
  "E8.0",
  "E8.25",
  "E8.5"
  # "mixed_gastrulation"
)


###############
## Load data ##
###############

# Load atlas metadata (which includes precomputed UMAP coordinates)
sample_metadata <- sample_metadata %>%
  .[stage%in%opts$stages] %>%
  .[,aggregated_celltype:=stringr::str_replace_all(celltype,opts$aggregated.celltypes)] %>%
  .[,aggregated_celltype:=factor(aggregated_celltype, levels=names(opts$celltype.colors))] %>%
  .[,stage:=factor(stage, levels=opts$stages)] %>%
  droplevels()

################
## Parse data ##
################

# Remove specific lineages
to.plot <- sample_metadata %>%
  # .[!celltype%in%c("ExE_ectoderm", "Parietal_endoderm","Visceral_endoderm","ExE_endoderm","PGC")] %>%
  # .[!celltype%in%c("ExE_ectoderm")] %>%
  .[,aggregated_celltype:=stringr::str_replace_all(aggregated_celltype,"_"," ")]
  
names(opts$celltype.colors) <- names(opts$celltype.colors) %>% stringr::str_replace_all("_"," ")
opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% unique(to.plot$aggregated_celltype)]

stopifnot(all(unique(to.plot$aggregated_celltype) %in% names(opts$celltype.colors)))
# unique(to.plot$aggregated_celltype)[!unique(to.plot$aggregated_celltype) %in% names(opts$celltype.colors)]

###################################
## Plot dimensionality reduction ##
###################################

# Colour by cell type
p <- ggplot(to.plot, aes(x=umapX, y=umapY)) +
  # geom_point(aes(colour=aggregated_celltype), size=0.1) +
  ggrastr::geom_point_rast(aes(colour=aggregated_celltype), size=0.05) +
  scale_color_manual(values=opts$celltype.colors) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

pdf(paste0(io$outdir,"/umap_per_celltype.pdf"), width=4.5, height=4.5)
print(p)
dev.off()


# Colour by stage
p <- ggplot(to.plot %>% setorder(stage), aes(x=umapX, y=umapY)) +
  # geom_point(aes(colour=stage), size=0.1) +
  ggrastr::geom_point_rast(aes(colour=stage), size=0.01) +
  scale_color_manual(values=opts$stage.colors) +
  guides(colour = guide_legend(override.aes = list(size=4.5))) +
  ggplot_theme_NoAxes() +
  theme(
    legend.title = element_blank(),
    legend.position = "right"
  )

pdf(paste0(io$outdir,"/umap_per_stage.pdf"), width=6.5, height=4.5)
print(p)
dev.off()

# Plot each sample separately
for (i in unique(to.plot$sample)) {
  print(i)
  
  to.plot_i <- to.plot %>% copy %>%
    .[,alpha:=1.0] %>%
    .[sample!=i,c("celltype","alpha"):=list("None",0.25)]
  
  p <- ggplot() +
    ggrastr::geom_point_rast(aes(x=umapX, y=umapY), size=1.5, color="grey", alpha=0.25, data=to.plot_i[sample!=i]) +
    ggrastr::geom_point_rast(aes(x=umapX, y=umapY, fill=aggregated_celltype), size=2, shape=21, alpha=1.0, data=to.plot_i[sample==i]) +
    scale_fill_manual(values=opts$celltype.colors) +
    theme_classic() +
    theme(
      legend.position="none"
    )
  
  pdf(sprintf("%s/umap_sample%s.pdf",io$outdir,i), width=5, height=5)
  print(p)
  dev.off()
}


#######################
## Plot trajectories ##
#######################

celltypes.to.highlight <- c("Haematoendothelial_progenitors","Blood_progenitors_1","Blood_progenitors_2","Erythroid1","Erythroid2","Erythroid3")
  
to.plot2 <- to.plot %>% copy %>%
  .[,alpha:=1.0] %>%
  .[!celltype%in%celltypes.to.highlight,c("celltype","alpha"):=list("None",0.25)]

# Colour by cell type
p <- ggplot(to.plot, aes(x=umapX, y=umapY)) +
  # geom_point(aes(colour=aggregated_celltype), size=0.1) +
  ggrastr::geom_point_rast(aes(colour=aggregated_celltype), size=0.05) +
  scale_color_manual(values=opts$celltype.colors) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  ggplot_theme_NoAxes() +
  theme(
    legend.position = "none"
  )

pdf(paste0(io$outdir,"/umap_per_celltype.pdf"), width=4.5, height=4.5)
print(p)
dev.off()