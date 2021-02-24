#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation10x/settings.R")
  io$mnn.script <- "/Users/ricard/gastrulation10x/mapping/mapping_mnn.R"
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/gastrulation10x/settings.R")
  io$mnn.script <- "/homes/ricard/gastrulation10x/mapping/mapping_mnn.R"
  io$tmpdir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/results/mapping/stages/tmp"
} 

# Cell types
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

# Stages
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

# Batches
opts$test_samples <- c(
  # they are all E7.5
  # "2",
  # "3",
  # "4",
  "6",
  "19",
  "20"
)

# Test mode (subset cells)?
opts$test <- FALSE


####################################################
## Run joint mapping with all cell types together ##
####################################################

opts$nPCs <- 50  # Number of PCs
opts$k <- 25     # Number of neighbours

for (i in opts$test_samples) {
  if (grepl("ricard",Sys.info()['nodename'])) {
    lsf <- ""
  } else {
    lsf <- sprintf("bsub -M 70000 -n 1 -q research-rh74 -o %s/%s.txt", io$tmpdir, i)
  }
  cmd <- sprintf("%s Rscript %s --test_samples %s --stages %s --nPCs %d --k %d", lsf, io$mnn.script, i, paste(opts$stages, collapse=" "), opts$nPCs, opts$k)
  if (isTRUE(opts$test)) cmd <- paste0(cmd, " --test")
  system(cmd)
}

##########################################
## Run mapping separately per cell type ##
##########################################

# opts$nPCs <- 5  # Number of PCs
# opts$k <- 10    # Number of neighbours
# 
# for (i in opts$test_samples) {
#   foo <- table(sample_metadata[sample==i,celltype])
#   celltypes <- names(foo[foo>50])
#   for (j in celltypes) {
#     if (grepl("ricard",Sys.info()['nodename'])) {
#       lsf <- ""
#     } else {
#       lsf <- sprintf("bsub -M 30000 -n 1 -q research-rh74 -o %s/%s_%s.txt", io$tmpdir, i,j)
#     }
#     cmd <- sprintf("%s Rscript %s --test_samples %s --celltypes %s --stages %s --nPCs %d --k %d", lsf, io$mnn.script, i, j, paste(opts$stages, collapse=" "), opts$nPCs, opts$k)
#     if (isTRUE(opts$test)) cmd <- paste0(cmd, " --test")
#     print(cmd); system(cmd)
#   }
# }
