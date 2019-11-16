library(Seurat)
srat <- readRDS("/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/processed/seurat/seurat.rds")
unique(srat@meta.data$stage)

# srat <- UpdateSeuratObject(srat)

# Subset E6.5 to E7.5 cells
cells <- rownames(srat@meta.data[srat@meta.data$stage %in% c("E6.5","E6.75","E7.0","E7.25","E7.5"),])
srat_subset <- subset(srat, cells=cells)
saveRDS(srat_subset, "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/processed/seurat_E6.5_to_E7.5.rds")

unique(srat_subset$meta.data$stage)

