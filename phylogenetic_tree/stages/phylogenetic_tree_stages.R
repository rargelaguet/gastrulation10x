library(ggpubr)

##############
## Settings ##
##############

source("/Users/ricard/gastrulation10x/settings.R")
io$outdir <- paste0(io$basedir,"/results/phylogenetic_tree/stages"); dir.create(io$outdir, showWarnings = F)

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

opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak"
  # "Caudal_epiblast",
  # "PGC",
  # "Anterior_Primitive_Streak",
  # "Notochord",
  # "Def._endoderm",
  # "Gut",
  # "Nascent_mesoderm",
  # "Mixed_mesoderm",
  # "Intermediate_mesoderm",
  # "Caudal_Mesoderm",
  # "Paraxial_mesoderm",
  # "Somitic_mesoderm",
  # "Pharyngeal_mesoderm",
  # "Cardiomyocytes",
  # "Allantois",
  # "ExE_mesoderm",
  # "Mesenchyme",
  # "Haematoendothelial_progenitors",
  # "Endothelium",
  # "Blood_progenitors_1",
  # "Blood_progenitors_2",
  # "Erythroid1",
  # "Erythroid2",
  # "Erythroid3",
  # "NMP",
  # "Rostral_neurectoderm",
  # "Caudal_neurectoderm",
  # "Neural_crest",
  # "Forebrain_Midbrain_Hindbrain",
  # "Spinal_cord",
  # "Surface_ectoderm",
  # "Visceral_endoderm",
  # "ExE_endoderm",
  # "ExE_ectoderm",
  # "Parietal_endoderm"
)

# Update metadata
sample_metadata <- sample_metadata %>% 
  .[stage%in%opts$stages & celltype%in%opts$celltypes]

###############
## Load data ##
###############

# Load SingleCellExperiment
sce <- readRDS(io$rna.sce)[,sample_metadata$cell]

# Filter genes
sce <- sce[rowSums(counts(sce))>25]

##################
## Computations ##
##################

# PCA
# pca <- prcomp(t(logcounts(sce)), rank.=5)
pca <- irlba::prcomp_irlba(t(logcounts(sce)), n=5)$x

dist.object <- dist(pca)

h <- hclust(dist.object, method = 'ward.D')
plot(h)


h2 <- hclust(matdistX, method = 'ward.D2')

plot(h2)


h2 <- hclust(matdistX, method = 'average')

plot(h2)

##########
## Save ##
##########

pdf(sprintf("%s/%s_affinity.pdf",io$outdir,i), width = 9, height = 5)
print(p)
dev.off()
fwrite(dt, paste0(io$outdir,"/celltype_affinity_E6.5_to_E7.5.txt.gz"))
