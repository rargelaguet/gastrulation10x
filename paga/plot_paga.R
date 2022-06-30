# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$outdir <- paste0(io$basedir,"/results/paga")


# Options
opts$celltypes = c(
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

#####################
## Load PAGA graph ##
#####################

source(here::here("paga/load_paga_graph.R"))

# Define cell type order
cellype.order <- rownames(connectivity.mtx)

##############################################
## Plot network, colour by cell type labels ##
##############################################

celltypes.to.highlight <- c("NMP","Caudal_Mesoderm","Somitic_mesoderm","Spinal_cord","Notochord","Primitive_Streak","Def._endoderm","Nascent_mesoderm")

p <- ggraph(igraph.paga.tbl, x = x, y = y) +
  # geom_edge_link0(aes(width = weight), edge_colour = "grey66", edge_alpha=0.2) +
  geom_edge_link(edge_colour = "grey66", edge_alpha=0.8, edge_width=0.25) +
  # geom_node_point(aes(fill=celltype), size=8, shape=21, alpha=1) +
  geom_node_point(aes(color=celltype), size=8, alpha=0.8) +
  # geom_node_text(aes(label=celltype), data=igraph.paga.tbl %>% filter(celltype%in%celltypes.to.highlight) %>% as_tibble) +
  # scale_fill_manual(values=opts$celltype.colors) +
  scale_color_manual(values=opts$celltype.colors) +
  theme_classic(base_size=14) +
  theme(
    axis.line = element_blank(), 
    axis.text = element_blank(),
    axis.ticks = element_blank(), 
    axis.title = element_blank(),
    legend.position = "none"
  )


pdf(paste0(io$outdir,"/paga_coloured_by_celltype.pdf"), width=4.25, height=4.5)
print(p)
dev.off()

##########################################
## Plot network, colour by trajectories ##
##########################################

trajectories <- list(
  "Mesoderm" = c("Epiblast", "Primitive_Streak", "Nascent_mesoderm"),
  "Endoderm" = c("Epiblast", "Primitive_Streak", "Def._endoderm", "Gut"),
  "Ectoderm" = c("Epiblast", "Rostral_neurectoderm","Forebrain_Midbrain_Hindbrain"),
  "Blood" = c("Haematoendothelial_progenitors", "Blood_progenitors_1", "Blood_progenitors_2", "Erythroid1", "Erythroid2", "Erythroid3")
)

i <- "Blood"
for (i in names(trajectories)) {
  
  opts$colors <- opts$celltype.colors[opts$celltypes]
  opts$colors[!names(opts$colors)%in%trajectories[[i]]] <- "gray70"
  
  opts$alpha = rep(1,length(opts$celltypes))
  names(opts$alpha) <- opts$celltypes
  opts$alpha[!names(opts$alpha)%in%trajectories[[i]]] <- 0.4
  
  opts$size = rep(6,length(opts$celltypes))
  names(opts$size) <- opts$celltypes
  opts$size[!names(opts$size)%in%trajectories[[i]]] <- 2
  
  opts$text_size = rep(2.3,length(opts$celltypes))
  names(opts$text_size) <- opts$celltypes
  opts$text_size[!names(opts$text_size)%in%trajectories[[i]]] <- 0
  
  # (...)
  
  pdf(sprintf("%s/trajectories/paga_coloured_by_%s_trajectory.pdf",io$outdir,i), width=5, height=4)
  print(p)
  dev.off()
  
}

##################################################
## Plot network, colour one cell type at a time ##
##################################################

for (i in opts$celltypes) {
  
  opts$colors <- opts$celltype.colors[opts$celltypes]
  # opts$colors[!names(opts$colors)==i] <- "gray70"
  
  opts$alpha = rep(1,length(opts$celltypes)); names(opts$alpha) <- opts$celltypes
  opts$alpha[!names(opts$alpha)==i] <- 0.30
  
  opts$size = rep(15,length(opts$celltypes)); names(opts$size) <- opts$celltypes
  opts$size[!names(opts$size)==i] <- 4.25
  
  opts$text_size = rep(4.25,length(opts$celltypes)); names(opts$text_size) <- opts$celltypes
  opts$text_size[!names(opts$text_size)==i] <- 1.75
  
  pdf(sprintf("%s/celltype/paga_coloured_by_celltype_%s.pdf",io$outdir,i), width=4.5, height=4)
  print(p)
  dev.off()
  
}