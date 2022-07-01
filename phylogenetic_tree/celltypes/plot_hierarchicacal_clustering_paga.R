##############
## Settings ##
##############

source("/Users/argelagr/gastrulation10x/settings.R")
source("/Users/argelagr/gastrulation10x/utils.R")

io$outdir <- paste0(io$basedir,"/results/phylogenetic_tree/celltypes"); dir.create(io$outdir, showWarnings = F)

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
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm"
  # "Visceral_endoderm",
  # "ExE_endoderm",
  # "ExE_ectoderm",
  # "Parietal_endoderm"
)

##############################
## Load connectivity matrix ##
##############################

paga_connectivity.mtx <- fread(io$paga.connectivity) %>% matrix.please

# Filter cell types
paga_connectivity.mtx <- paga_connectivity.mtx[opts$celltypes,opts$celltypes]

# Parse
paga_distance.mtx  <- 1 - paga_connectivity.mtx
diag(paga_distance.mtx) <- 0

#############################
## Hierarchical clustering ##
#############################

paga_distance.dist <- as.dist(paga_distance.mtx)

h <- hclust(paga_distance.dist, method = 'ward.D')
plot(h)


h2 <- hclust(paga_distance.dist, method = 'ward.D2')
plot(h2)


h2 <- hclust(paga_distance.dist, method = 'average')
plot(h2)
