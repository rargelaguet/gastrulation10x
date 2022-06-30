##############################################
## Plot network, colour by cell type labels ##
##############################################

p <- ggnet2(
  net = net,
  mode = c("x", "y"),
  color = opts$celltype.colors[opts$celltypes],
  node.size = 6,
  edge.size = 0.15,
  edge.color = "grey",
  label = TRUE,
  label.size = 2.3
)

pdf(paste0(io$outdir,"/paga_coloured_by_celltype.pdf"), width=4.5, height=4)
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
  
  p <- ggnet2(
    net = net,
    mode = c("x", "y"),
    color = opts$colors,
    # color = opts$celltype.colors[opts$celltypes],
    node.size = opts$size,
    edge.size = 0.10,
    edge.color = "grey",
    alpha = opts$alpha,
    label = TRUE,
    label.size = opts$text_size
  ) + guides(size="none")
  
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
  
  p <- ggnet2(
    net = net,
    mode = c("x", "y"),
    color = opts$colors,
    # color = opts$celltype.colors[opts$celltypes],
    node.size = opts$size,
    edge.size = 0.10,
    edge.color = "grey",
    alpha = opts$alpha,
    label = TRUE,
    label.size = opts$text_size
  ) + guides(size="none")
  
  pdf(sprintf("%s/celltype/paga_coloured_by_celltype_%s.pdf",io$outdir,i), width=4.5, height=4)
  print(p)
  dev.off()
  
}


########################################################
## Plot network, colour multiple cell types at a time ##
########################################################

celltypes <- list(
  "blood_progenitors" = c("Blood_progenitors_1", "Blood_progenitors_2"),
  "erythroids" = c("Erythroid1", "Erythroid2", "Erythroid3")
)

for (i in names(celltypes)) {
  
  opts$colors <- opts$celltype.colors[opts$celltypes]
  # opts$colors[!names(opts$colors)==i] <- "gray70"
  
  opts$alpha = rep(1,length(opts$celltypes)); names(opts$alpha) <- opts$celltypes
  opts$alpha[!names(opts$alpha)%in%celltypes[[i]]] <- 0.30
  
  opts$size = rep(15,length(opts$celltypes)); names(opts$size) <- opts$celltypes
  opts$size[!names(opts$size)%in%celltypes[[i]]] <- 4.25
  
  opts$text_size = rep(4.25,length(opts$celltypes)); names(opts$text_size) <- opts$celltypes
  opts$text_size[!names(opts$text_size)%in%celltypes[[i]]] <- 1.75
  
  p <- ggnet2(
    net = net,
    mode = c("x", "y"),
    color = opts$colors,
    # color = opts$celltype.colors[opts$celltypes],
    node.size = opts$size,
    edge.size = 0.10,
    edge.color = "grey",
    alpha = opts$alpha,
    label = TRUE,
    label.size = opts$text_size
  ) + guides(size="none")
  
  pdf(sprintf("%s/celltype/paga_coloured_by_celltype_%s.pdf",io$outdir,i), width=4.5, height=4)
  print(p)
  dev.off()
  
}
