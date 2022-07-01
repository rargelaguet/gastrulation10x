# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

io$trajectory.dir <- file.path(io$basedir,"results/trajectories/blood_scanpy")
io$sce <- file.path(io$basedir,"results/trajectories/blood_scanpy/blood_SingleCellExperiment.rds")
io$outdir <- file.path(io$basedir,"results/trajectories/blood_scanpy/stuff")
dir.create(io$outdir, showWarnings=F)

##########################
## Load sample metadata ##
##########################

sample_metadata.dt <- fread(file.path(io$trajectory.dir,"blood_sample_metadata.txt.gz")) %>%
  .[,c("cell","stage","celltype")]

#####################
## Load trajectory ##
#####################

trajectory.dt <- fread(file.path(io$trajectory.dir,"blood_trajectory.txt.gz")) %>%
  setnames(c("FA1","FA2"),c("V1","V2"))

#########################
## Load RNA expression ##
#########################

sce <- load_SingleCellExperiment(io$sce, normalise = T, cells=sample_metadata.dt$cell)

# Add sample metadata as colData
colData(sce) <- sample_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

# Rename genes
# gene_metadata <- fread(io$gene_metadata) %>% .[ens_id%in%rownames(sce)]
# foo <- gene_metadata$symbol
# names(foo) <- gene_metadata$ens_id
# sce <- sce[rownames(sce) %in% names(foo),]
# rownames(sce) <- foo[rownames(sce)]

# Save SingleCellExperiment
# assays(sce)[["logcounts"]] <- NULL
# saveRDS(sce, file.path(io$outdir,"blood_SingleCellExperiment.rds"))

################
## Parse data ##
################

# names(opts$celltype.colors) <- names(opts$celltype.colors) %>% stringr::str_replace_all("_"," ")
# opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% unique(to.plot$aggregated_celltype)]

###################################
## Plot dimensionality reduction ##
###################################

to.plot <- trajectory.dt %>% merge(sample_metadata.dt, by="cell")

to.plot.subset <- to.plot[,.SD[sample.int(1000)],by="celltype"]

p <- ggplot(to.plot.subset, aes(x=V1, y=V2, fill=celltype)) +
  geom_point(size=1, shape=21, stroke=0.1) +
  ggrastr::geom_point_rast(size=1, shape=21, stroke=0.1) +
  scale_fill_manual(values=opts$celltype.colors[names(opts$celltype.colors)%in%unique(to.plot$celltype)]) +
  # ggrepel::geom_text_repel(aes_string(label="celltype"), data=to.plot[,.(V1=median(V1), V2=median(V2)), by="celltype"], size=4) +
  theme_classic() +
  labs(x="", y="") +
  guides(fill = guide_legend(override.aes = list(size=3))) +
  ggplot_theme_NoAxes() +
  theme(
    legend.position="right"
  )

pdf(paste0(io$outdir,"/blood_trajectory_forcelayout_per_celltype_legend.pdf"), width=4, height=4)
# pdf(paste0(io$outdir,"/legend.pdf"), width=8, height=4)
print(p)
dev.off()


##########################
## Plot gene expression ##
##########################

# Select genes
# genes.to.plot <- c("Dnmt1","Dnmt3a","Dnmt3b","Tet1","Tet2","Tet3")
genes.to.plot <- c("Hbb-y","Hbb-bh1","Hba-x","Hba-a1","Hba-a2")
# genes.to.plot <- grep("Hb",rownames(sce),value=T)
sce_filt <- sce[genes.to.plot]

# Subset cells
sce_filt <- sce_filt[,sample.int(n=ncol(sce_filt), size=ncol(sce_filt)/1.5)]

# Denoise
pca.mtx <- trajectory.dt %>% matrix.please %>% .[colnames(sce_filt),]
rna_denoised.mtx <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(sce_filt)), D=pdist(pca.mtx), k=25)
colnames(rna_denoised.mtx) <- colnames(sce_filt)
assay(sce_filt,"logcounts_denoised") <- rna_denoised.mtx
rm(rna_denoised.mtx); gc(reset=T)

# Plot
for (i in genes.to.plot) {
  
  to.plot <- data.table(
    cell = colnames(sce_filt),
    expr = logcounts(sce_filt)[i,],
    expr_denoised = assay(sce_filt,"logcounts_denoised")[i,]
  ) %>% 
    melt(id.vars=c("cell"), value.name="expr") %>%
    merge(trajectory.dt, by="cell")
  
  p <- ggplot(to.plot[variable=="expr_denoised"], aes(x=V1, y=V2)) +
    geom_point(color="gray90", size=0.75, alpha=0.5, data=to.plot[expr==0]) +
    geom_point(aes(fill=expr), size=1, shape=21, stroke=0.02, data=to.plot[expr>0]) +
    scale_fill_gradient2(low = "gray50", mid="gray90", high = "darkred") +
    # facet_wrap(~variable, nrow=1) +
    # ggrepel::geom_text_repel(aes_string(label="celltype"), data=to.plot[,.(V1=median(V1), V2=median(V2)), by="celltype"]) +
    theme_classic() +
    labs(x="", y="", title=i) +
    ggplot_theme_NoAxes() +
    theme(
      legend.position="none"
    )
  
  pdf(sprintf("%s/blood_trajectory_expr_%s.pdf",io$outdir,i), width=8, height=6.5)
  print(p)
  dev.off()
    
}
