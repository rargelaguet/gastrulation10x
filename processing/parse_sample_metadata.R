source("/Users/ricard/gastrulation10x/settings.R")

cols.to.select <- c("cell", "barcode", "sample", "stage", "sequencing.batch", "doub.density", "doublet", "stripped", "celltype", "umapX", "umapY")
sample_metadata <- sample_metadata[,cols.to.select,with=F]

to.merge <- c(
  "Caudal neurectoderm" = "Neurectoderm",
  "Rostral neurectoderm" = "Neurectoderm",
  "Erythroid1" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid3" = "Erythroid",
  "Blood progenitors 1" = "Blood progenitors",
  "Blood progenitors 2" = "Blood progenitors",
  "Intermediate mesoderm" = "Mixed mesoderm"
)

sample_metadata <- sample_metadata %>%
  .[,celltype2:=celltype] %>%
  .[,celltype2:=stringr::str_replace_all(celltype2, to.merge)]

fwrite(sample_metadata, io$metadata, sep="\t", na="NA", quote=F)

