source("/Users/ricard/gastrulation10x/settings.R")

# cols.to.select <- c("cell", "barcode", "sample", "stage", "sequencing.batch", "doub.density", "doublet", "stripped", "celltype", "celltype2", "umapX", "umapY")
# sample_metadata <- sample_metadata[,cols.to.select,with=F]

# to.merge <- c(
#   "Caudal_neurectoderm" = "Neurectoderm",
#   "Rostral_neurectoderm" = "Neurectoderm",
#   "Erythroid1" = "Erythroid",
#   "Erythroid2" = "Erythroid",
#   "Erythroid3" = "Erythroid",
#   "Blood_progenitors_1" = "Blood_progenitors",
#   "Blood_progenitors_2" = "Blood_progenitors",
#   "Intermediate_mesoderm" = "Mixed_mesoderm"
# )

to.merge <- c(
  "Epiblast" = "Epiblast-PS",
  "Primitive_Streak" = "Epiblast-PS",
  "ExE_ectoderm" = "ExE_ectoderm",
  "Visceral_endoderm" = "ExE_endoderm",
  "ExE_endoderm" = "ExE_endoderm",
  "Nascent_mesoderm" = "Nascent_mesoderm",
  "Neurectoderm" = "Neuroectoderm",
  "Blood_progenitors" = "Blood_progenitors",
  "Mixed_mesoderm" = "Mesoderm",
  "ExE_mesoderm" = "ExE_mesoderm",
  "Pharyngeal_mesoderm" = "Mesoderm",
  "Caudal_epiblast" = "Caudal_epiblast",
  "PGC" = "PGC",
  "Mesenchyme" = "Mesenchyme",
  "Haematoendothelial_progenitors" = "Blood_progenitors",
  "Surface_ectoderm" = "Surface_ectoderm",
  "Gut" = "Endoderm",
  "Paraxial_mesoderm" = "Mesoderm",
  "Notochord" = "Notochord",
  "Somitic_mesoderm" = "Mesoderm",
  "Caudal_Mesoderm" = "Mesoderm",
  "Erythroid" = "Erythroid",
  "Def._endoderm" = "Endoderm",
  "Parietal_endoderm" = "Parietal_endoderm",
  "Allantois" = "ExE_mesoderm",
  "Anterior_Primitive_Streak" = "Epiblast-PS",
  "Endothelium" = "Endothelium",
  "Forebrain_Midbrain_Hindbrain" = "Neuroectoderm",
  "Spinal_cord" = "Spinal_cord",
  "Cardiomyocytes" = "Cardiomyocytes",
  "NMP" = "NMP",
  "Neural_crest" = "Neuroectoderm"
)

sample_metadata <- sample_metadata %>%
  .[,celltype3:=celltype2] %>%
  .[,celltype3:=stringr::str_replace_all(celltype3, to.merge)]

table(sample_metadata$celltype3)
fwrite(sample_metadata, io$metadata, sep="\t", na="NA", quote=F)

