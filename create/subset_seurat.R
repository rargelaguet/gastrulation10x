library(Seurat)
srat <- readRDS("/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/processed/seurat.rds")

# Subset E6.5 cells
# cells <- rownames(srat@meta.data[srat@meta.data$stage == "E6.5",])
# srat_subset <- SubsetData(srat, cells.use=cells, do.center=T, do.scale=T, do.clean=T)
# saveRDS(srat_subset, "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/atlas/processed/seurat_E6.5.rds")

# Subset E7.0 cells
# cells <- rownames(srat@meta.data[srat@meta.data$stage == "E7.0",])
# srat_subset <- SubsetData(srat, cells.use=cells, do.center=T, do.scale=T, do.clean=T)
# saveRDS(srat_subset, "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/atlas/processed/seurat_E7.0.rds")

# Subset E7.5 cells
# cells <- rownames(srat@meta.data[srat@meta.data$stage == "E7.5",])
# srat_subset <- SubsetData(srat, cells.use=cells, do.center=T, do.scale=T, do.clean=T)
# saveRDS(srat_subset, "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/atlas/processed/seurat_E7.5.rds")


# Subset E6.5 to E7.5 cells
cells <- rownames(srat@meta.data[srat@meta.data$stage %in% c("E6.5","E7.0","E7.5"),])
srat_subset <- SubsetData(srat, cells.use=cells, do.center=T, do.scale=T, do.clean=T)
saveRDS(srat_subset, "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/processed/seurat_E6.5_to_E7.5.rds")


