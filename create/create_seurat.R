library(data.table)
library(purrr)
# library(scater)
library(Seurat)

directory <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/atlas"

# Load cell metadata
barcode.loc <- file.path(directory, "meta.tab")
cell.info <- read.delim(barcode.loc, header = TRUE, colClasses = "character", stringsAsFactors = FALSE, sep="\t")
rownames(cell.info) <- cell.info$cell
  
# Load gene metadata
gene.loc <- file.path(directory, "genes.tsv")
gene.info <- read.delim(gene.loc, header = FALSE, colClasses = "character", stringsAsFactors = FALSE, sep="\t")
colnames(gene.info) <- c("ens_id","gene_id")

# Load matrix  
matrix.loc <- file.path(directory, "raw_counts.mtx")
data_mat <- Matrix::readMM(matrix.loc)
rownames(data_mat) <- gene.info$gene_id
colnames(data_mat) <- cell.info$cell

# Subset  
data_mat <- data_mat[,cell.info$doublet==F & cell.info$stripped==F]
cell.info <- cell.info[cell.info$doublet==F & cell.info$stripped==F,]

# Create seurat object
srat <- CreateSeuratObject(
  raw.data = data_mat,
  min.genes = 1, min.cells = 500,
  meta.data = cell.info,
  do.scale=FALSE, do.center=FALSE
)
rm(data_mat)

# Normalize data
srat <- NormalizeData(srat, scale.factor=1000)

# Scale data
srat <- ScaleData(srat, 
  vars.to.regress="nUMI",
  model.use = "linear",
  do.scale = FALSE,
  do.center = TRUE, 
  do.par = TRUE, num.cores = 4
)

# Save
saveRDS(srat, "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/atlas/seurat.rds")