celltype_colours = c("Epiblast" = "#635547",
                     "Primitive Streak" = "#DABE99",
                     "Caudal epiblast" = "#9e6762",
                     
                     "PGC" = "#FACB12",
                     
                     "Anterior Primitive Streak" = "#c19f70",
                     "Notochord" = "#0F4A9C",
                     "Def. endoderm" = "#F397C0",
                     "Gut" = "#EF5A9D",
                     
                     "Nascent mesoderm" = "#C594BF",
                     "Mixed mesoderm" = "#DFCDE4",
                     "Intermediate mesoderm" = "#139992",
                     "Caudal Mesoderm" = "#3F84AA",
                     "Paraxial mesoderm" = "#8DB5CE",
                     "Somitic mesoderm" = "#005579",
                     "Pharyngeal mesoderm" = "#C9EBFB",
                     "Cardiomyocytes" = "#B51D8D",
                     "Allantois" = "#532C8A",
                     "ExE mesoderm" = "#8870ad",
                     "Mesenchyme" = "#cc7818",
                     
                     "Haematoendothelial progenitors" = "#FBBE92",
                     "Endothelium" = "#ff891c",
                     "Blood progenitors 1" = "#f9decf",
                     "Blood progenitors 2" = "#c9a997",
                     "Erythroid1" = "#C72228",
                     "Erythroid2" = "#f79083",
                     "Erythroid3" = "#EF4E22",
                     
                     "NMP" = "#8EC792",
                     
                     "Rostral neurectoderm" = "#65A83E",
                     "Caudal neurectoderm" = "#354E23",
                     "Neural crest" = "#C3C388",
                     "Forebrain/Midbrain/Hindbrain" = "#647a4f",
                     "Spinal cord" = "#CDE088",
                     
                     "Surface ectoderm" = "#f7f79e",
                     
                     "Visceral endoderm" = "#F6BFCB",
                     "ExE endoderm" = "#7F6874",
                     "ExE ectoderm" = "#989898",
                     "Parietal endoderm" = "#1A1A1A"
                     
)


srat <- readRDS("/Users/ricard/data/gastrulation10x/atlas/processed/seurat_E6.5.rds")
meta.data <- srat@meta.data

# Remove genes
srat <- CreateSeuratObject( srat@raw.data[apply(srat@data,1,var)>0,] )
srat <- CreateSeuratObject( srat@raw.data[rownames(srat@raw.data) != rownames(srat@data)[duplicated(rownames(srat@data))],] )
srat <- AddMetaData(srat, metadata = meta.data)

srat <- NormalizeData(srat, scale.factor=1000)

srat <- ScaleData(
  srat, vars.to.regress=c("sample"),
  model.use = "linear",
  do.scale = FALSE, do.center = TRUE, 
  do.par = TRUE, num.cores = 4
)

srat <- FindVariableGenes(srat, x.low.cutoff = 0.01, x.high.cutoff = 1.0)

# Run
srat <- RunPCA(srat, pcs.compute = 15)
# srat <- RunTSNE(srat) # run t-SNE after pre-computing PCA
srat <- RunUMAP(srat) # run t-SNE after pre-computing PCA

# Plot
# DimPlot(srat, reduction.use="pca", group.by = "celltype", dim.1 = 2, dim.2 = 3)
# DimPlot(srat, reduction.use="tsne", group.by = "celltype")
DimPlot(srat, reduction.use="umap", group.by = "celltype", cols.use=celltype_colours)


