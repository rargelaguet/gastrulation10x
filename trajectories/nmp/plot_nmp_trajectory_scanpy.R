# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

library(scran)
library(scater)

#####################
## Define settings ##
#####################

io$trajectory.dir <- file.path(io$basedir,"results/trajectories/nmp_somitic_spinal")
io$outdir <- file.path(io$basedir,"results/trajectories/nmp_somitic_spinal"); dir.create(io$outdir, showWarnings=F)

##########################
## Load sample metadata ##
##########################

sample_metadata.dt <- fread(file.path(io$trajectory.dir,"nmp_sample_metadata.txt.gz")) %>%
  .[,c("cell","stage","celltype")]

#####################
## Load trajectory ##
#####################

trajectory.dt <- fread(file.path(io$trajectory.dir,"nmp_trajectory.txt.gz")) %>%
  setnames(c("FA1","FA2"),c("V1","V2"))

#########################
## Load RNA expression ##
#########################

sce <- load_SingleCellExperiment(io$sce, normalise = T, cells=sample_metadata.dt$cell)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

# Rename genes
gene_metadata <- fread(io$gene_metadata) %>% .[ens_id%in%rownames(sce)]
foo <- gene_metadata$symbol
names(foo) <- gene_metadata$ens_id
sce <- sce[rownames(sce) %in% names(foo),]
rownames(sce) <- foo[rownames(sce)]

###################################
## Plot dimensionality reduction ##
###################################

to.plot <- trajectory.dt %>% merge(sample_metadata.dt, by="cell")

p <- ggplot(to.plot, aes(x=V1, y=V2, fill=celltype)) +
  # geom_point(size=1, shape=21, stroke=0.2) +
  ggrastr::geom_point_rast(size=1, shape=21, stroke=0.2, raster.dpi=150) +
  scale_fill_manual(values=opts$celltype.colors) +
  # ggrepel::geom_text_repel(aes_string(label="celltype"), data=to.plot[,.(V1=median(V1), V2=median(V2)), by="celltype"]) +
  theme_classic() +
  labs(x="", y="") +
  ggplot_theme_NoAxes() +
  theme(
    legend.position="none"
  )

pdf(paste0(io$outdir,"/nmp_trajectory_color_per_celltype.pdf"), width=5, height=5)
print(p)
dev.off()

#############
## Denoise ##
#############

# Feature selection
sce <- sce[rowSums(counts(sce))>=50,]
decomp <- modelGeneVar(sce)
hvgs <- decomp[order(decomp$FDR),] %>% head(n=500) %>% rownames

# dimensionality reduction
sce.filt <- runPCA(sce[hvgs,], ncomponents = 25)

# kNN denoise
rna_denoised.mtx <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(sce)), D=pdist(reducedDim(sce.filt,"PCA")), k=15)
colnames(rna_denoised.mtx) <- colnames(sce.filt)
assay(sce,"logcounts_denoised") <- rna_denoised.mtx
rm(rna_denoised.mtx); gc(reset=T)

##########################
## Plot gene expression ##
##########################

genes.to.plot <- c("Sox2","T","Hoxb9")

for (i in genes.to.plot) {
  
  to.plot <- data.table(
    cell = colnames(sce),
    # expr = logcounts(sce)[i,]
    expr = assay(sce,"logcounts_denoised")[i,]
  ) %>% merge(trajectory.dt, by="cell")
  
  p <- ggplot(to.plot, aes(x=V1, y=V2, fill=expr)) +
    geom_point(size=1, shape=21, stroke=0.2) +
    scale_fill_gradient2(low = "gray50", mid="gray90", high = "red") +
    # ggrepel::geom_text_repel(aes_string(label="celltype"), data=to.plot[,.(V1=median(V1), V2=median(V2)), by="celltype"]) +
    theme_classic() +
    labs(x="", y="", title=i) +
    ggplot_theme_NoAxes() +
    theme(
      legend.position="none"
    )
  
  pdf(sprintf("%s/nmp_trajectory_color_by_%s.pdf",io$outdir,i), width=7, height=5)
  print(p)
  dev.off()
    
}

#########################
## Plot T/Sox2 overlap ##
#########################

to.plot <- data.table(
  cell = colnames(sce),
  # expr = logcounts(sce)[i,]
  expr = minmax.normalisation(assay(sce,"logcounts_denoised")["Sox2",]*assay(sce,"logcounts_denoised")["T",])
) %>% merge(trajectory.dt, by="cell")

p <- ggplot(to.plot, aes(x=V1, y=V2, fill=expr)) +
  geom_point(size=1, shape=21, stroke=0.2) +
  scale_fill_gradient2(low = "gray50", mid="gray90", high = "red") +
  # ggrepel::geom_text_repel(aes_string(label="celltype"), data=to.plot[,.(V1=median(V1), V2=median(V2)), by="celltype"]) +
  theme_classic() +
  labs(x="", y="", title="T*Sox2") +
  ggplot_theme_NoAxes() +
  theme(
    legend.position="none"
  )

pdf(sprintf("%s/nmp_trajectory_color_by_T_Sox2_product.pdf",io$outdir), width=7, height=5)
print(p)
dev.off()

##########
## Save ##
##########

saveRDS(sce, file.path(io$outdir,"nmp_SingleCellExperiment.rds"))
