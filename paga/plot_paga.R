library(GGally)
library(network)
library(sna)

#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation10x/settings.R")
  source("/Users/ricard/gastrulation10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation10x/settings.R")
  source("/homes/ricard/gastrulation10x/utils.R")
} else {
  stop("Computer not recognised")
}

# I/O
io$outdir <- paste0(io$basedir,"/results/paga")


# Options
opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  # "PGC",
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
  # "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)

#######################
## Load network data ##
#######################

df.connectivity <- fread(io$paga.connectivity) %>%
  matrix.please %>% .[opts$celltypes,opts$celltypes]

df.coordinates <- fread(io$paga.coordinates) %>% 
  matrix.please %>% .[opts$celltypes,]

####################
## Create network ##
####################

# Parse data
df.connectivity[df.connectivity<0.20] <- 0

# Create network
net = network(df.connectivity)

# Define coordinates
x = gplot.layout.fruchtermanreingold(net, NULL)
net %v% "x" = df.coordinates[, 1]
net %v% "y" = df.coordinates[, 2]

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

pdf(paste0(io$outdir,"/paga_coloured_by_celltype.pdf"), width=6, height=4)
print(p)
dev.off()


##########################################
## Plot network, colour by trajectories ##
##########################################

trajectories <- list(
  "Mesoderm" = c("Epiblast", "Primitive_Streak", "Nascent_mesoderm", "Mixed_mesoderm", "Pharyngeal_mesoderm"),
  "Endoderm" = c("Epiblast", "Primitive_Streak", "Anterior_Primitive_Streak", "Def._endoderm", "Gut"),
  "Ectoderm" = c("Epiblast", "Rostral_neurectoderm","Forebrain_Midbrain_Hindbrain"),
  "Blood" = c("Haematoendothelial_progenitors", "Blood_progenitors_1", "Blood_progenitors_2", "Erythroid1", "Erythroid2", "Erythroid3")
)

for (i in names(trajectories)) {
  
  opts$colors <- opts$celltype.colors[opts$celltypes]
  opts$colors[!names(opts$colors)%in%trajectories[[i]]] <- "gray70"
  
  opts$alpha = rep(1,length(opts$celltypes))
  names(opts$alpha) <- opts$celltypes
  opts$alpha[!names(opts$alpha)%in%trajectories[[i]]] <- 0.4
  
  opts$size = rep(6,length(opts$celltypes))
  names(opts$size) <- opts$celltypes
  opts$size[!names(opts$size)%in%trajectories[[i]]] <- 3
  
  opts$text_size = rep(2.3,length(opts$celltypes))
  names(opts$text_size) <- opts$celltypes
  opts$text_size[!names(opts$text_size)%in%trajectories[[i]]] <- 1.5
  
  p <- ggnet2(
    net = net,
    mode = c("x", "y"),
    color = opts$colors,
    # color = opts$celltype.colors[opts$celltypes],
    node.size = opts$size,
    edge.size = 0.15,
    edge.color = "grey",
    alpha = opts$alpha,
    label = TRUE,
    label.size = opts$text_size
  ) + guides(size=F)
  
  pdf(sprintf("%s/trajectories/paga_coloured_by_%s_trajectory.pdf",io$outdir,i), width=5, height=4)
  print(p)
  dev.off()
  
}
