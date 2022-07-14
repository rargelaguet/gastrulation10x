here::i_am("dimensionality_reduction/precomputed_dimensionality_reduction.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

io$outdir <- paste0(io$basedir,"/results/dimensionality_reduction/pdf")

opts$celltypes <- c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)

opts$rename_celltypes <- c(
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

opts$subset_cells <- FALSE

###############
## Load data ##
###############

# Load atlas metadata (which includes precomputed UMAP coordinates)
sample_metadata <- fread(io$metadata) %>%
  .[stripped==F & doublet==F] %>%
  .[stage%in%opts$stages & celltype%in%opts$celltypes] #%>%
  # .[,aggregated_celltype:=stringr::str_replace_all(celltype,opts$aggregated.celltypes)] %>%
  # .[,aggregated_celltype:=factor(aggregated_celltype, levels=names(opts$celltype.colors))] %>%
  # .[,stage:=factor(stage, levels=opts$stages)] %>%
  # droplevels()

if (opts$subset_cells) sample_metadata <- sample_metadata[sample(.N,7e4)]

################
## Parse data ##
################

to.plot <- sample_metadata %>%
  # .[,celltype:=stringr::str_replace_all(celltype,opts$rename_celltypes)] %>%
  # .[,celltype:=factor(celltype, levels=names(opts$celltype.colors))] %>%
  # .[!celltype%in%c("ExE_ectoderm", "Parietal_endoderm","Visceral_endoderm","ExE_endoderm","PGC")] %>%
  # .[!celltype%in%c("ExE_ectoderm")] %>%
  .[,celltype:=stringr::str_replace_all(celltype,"_"," ")]
  
names(opts$celltype.colors) <- names(opts$celltype.colors) %>% stringr::str_replace_all("_"," ")
opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% unique(to.plot$celltype)]

stopifnot(all(unique(to.plot$celltype) %in% names(opts$celltype.colors)))
# unique(to.plot$celltype)[!unique(to.plot$celltype) %in% names(opts$celltype.colors)]

#########################################################
## Plot dimensionality reduction coloured by cell type ##
#########################################################

# to.plot <- to.plot[,.SD[sample.int(1000)],by="celltype"]
# to.plot.subset <- to.plot[,.SD[sample.int(5e4)]]
# to.plot.subset <- to.plot

p <- ggplot(to.plot, aes(x=umapX, y=umapY)) +
  # geom_point(aes(colour=celltype), size=0.1) +
  ggrastr::geom_point_rast(size=0.05, aes(colour=celltype), alpha=0.80, raster.dpi=150) +
  ggrepel::geom_text_repel(aes_string(label="celltype"), size=4, data=to.plot[,.(umapX=median(umapX), umapY=median(umapY)), by="celltype"]) +
  scale_color_manual(values=opts$celltype.colors) +
  guides(colour = guide_legend(nrow=4, byrow=TRUE, override.aes = list(size=3.5))) +
  theme_classic() +
  ggplot_theme_NoAxes() +
  theme(
    legend.position = "none"
  )

pdf(paste0(io$outdir,"/umap_per_celltype.pdf"), width=7, height=6.25)
print(p)
dev.off()


#####################################################
## Plot dimensionality reduction coloured by stage ##
#####################################################

p <- ggplot(to.plot %>% setorder(-stage), aes(x=umapX, y=umapY)) +
  # geom_point(aes(colour=stage), size=0.1) +
  ggrastr::geom_point_rast(aes(colour=stage), size=0.05, alpha=0.80, raster.dpi=150) +
  scale_color_manual(values=opts$stage.colors) +
  guides(colour = guide_legend(override.aes = list(size=4.5))) +
  ggplot_theme_NoAxes() +
  theme(
    legend.title = element_blank(),
    legend.position = "none"
  )

pdf(paste0(io$outdir,"/umap_per_stage.pdf"), width=7, height=6.1)
print(p)
dev.off()

###############################
## Plot one sample at a time ##
###############################

# for (i in unique(to.plot$sample)) {
#   print(i)
#   
#   to.plot_i <- to.plot %>% copy %>%
#     .[,alpha:=1.0] %>%
#     .[sample!=i,c("celltype","alpha"):=list("None",0.25)]
#   
#   p <- ggplot() +
#     ggrastr::geom_point_rast(aes(x=umapX, y=umapY), size=1.5, color="grey", alpha=0.25, data=to.plot_i[sample!=i]) +
#     ggrastr::geom_point_rast(aes(x=umapX, y=umapY, fill=celltype), size=2, shape=21, alpha=1.0, data=to.plot_i[sample==i]) +
#     scale_fill_manual(values=opts$celltype.colors) +
#     theme_classic() +
#     ggplot_theme_NoAxes() +
#     theme(
#       legend.position="none"
#     )
#   
#   pdf(sprintf("%s/umap_sample%s.pdf",io$outdir,i), width=5, height=5)
#   print(p)
#   dev.off()
# }


#################################
## Plot one celltype at a time ##
#################################

celltypes.to.plot <- unique(to.plot$celltype)

for (i in celltypes.to.plot) {
  print(i)

  to.plot_i <- to.plot %>% copy %>%
    .[,alpha:=1.0] %>%
    .[!celltype%in%i,c("celltype","alpha"):=list("None",0.25)]

  p <- ggplot() +
    ggrastr::geom_point_rast(aes(x=umapX, y=umapY), size=0.25, color="grey", alpha=0.25, raster.dpi=150, data=to.plot_i[!celltype%in%i][sample.int(.N,5e4)]) +
    ggrastr::geom_point_rast(aes(x=umapX, y=umapY, fill=celltype), size=0.80, shape=21, alpha=1.0, stroke=0.075, raster.dpi=150, data=to.plot_i[celltype%in%i]) +
    scale_fill_manual(values=opts$celltype.colors) +
    theme_classic() +
    ggplot_theme_NoAxes() +
    theme(
      legend.position="none"
    )

  pdf(sprintf("%s/umap_coloured_by_%s.pdf",io$outdir,gsub(" ","-",i)), width=6, height=6)
  print(p)
  dev.off()
}

#######################
## Plot trajectories ##
#######################

# celltypes.to.highlight <- c("Haematoendothelial progenitors","Blood progenitors 1","Blood progenitors 2","Erythroid1","Erythroid2","Erythroid3")
# celltypes.to.highlight <- c("Erythroid1","Erythroid2","Erythroid3")
celltypes.to.highlight <- c("Blood progenitors 1","Blood progenitors 2")
# celltypes.to.highlight <- c("ExE ectoderm","ExE endoderm","Parietal endoderm")
  
to.plot2 <- to.plot %>% copy %>%
  .[,alpha:=1.0] %>%
  .[!celltype%in%celltypes.to.highlight,c("celltype","alpha"):=list("None",0.25)]

# Colour by cell type
p <- ggplot(to.plot2, aes(x=umapX, y=umapY)) +
  ggrastr::geom_point_rast(aes(x=umapX, y=umapY), size=0.25, color="grey", alpha=0.25, raster.dpi=150, data=to.plot2[!celltype%in%celltypes.to.highlight][sample.int(.N,5e4)]) +
  ggrastr::geom_point_rast(aes(x=umapX, y=umapY, fill=celltype), size=0.80, shape=21, alpha=1.0, stroke=0.1, raster.dpi=150, data=to.plot2[celltype%in%celltypes.to.highlight]) +
  scale_fill_manual(values=opts$celltype.colors[celltypes.to.highlight]) +
  ggplot_theme_NoAxes() +
  theme(
    legend.position = "none"
  )

pdf(paste0(io$outdir,"/umap_coloured_by_blood_progenitors.pdf"), width=6, height=6)
print(p)
dev.off()

##########
## Save ##
##########

umap.dt <- sample_metadata %>%
  .[,c("cell","umapX","umapY","celltype")] %>%
  setnames(c("umapX","umapY"),c("V1","V2")) %>%
  .[,c("V1","V2"):=list(round(V1,2),round(V2,2))]
fwrite(umap.dt, file.path(io$outdir,"umap_coordinates.txt.gz"), sep="\t", quote=F, na="NA")
