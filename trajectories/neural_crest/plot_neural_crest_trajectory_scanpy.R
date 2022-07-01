# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

io$trajectory.dir <- file.path(io$basedir,"results/trajectories/blood_scanpy")
io$outdir <- file.path(io$basedir,"results/trajectories/blood_scanpy/stuff")
dir.create(io$outdir, showWarnings=F)

##########################
## Load sample metadata ##
##########################

sample_metadata.dt <- fread(file.path(io$trajectory.dir,"sample_metadata.txt.gz")) %>%
  .[,c("cell","stage","celltype")]

########################
## Load gene metadata ##
########################

gene_metadata <- fread(io$gene_metadata) %>%
  .[ens_id%in%rownames(sce)]

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
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

# Rename genes
foo <- gene_metadata$symbol
names(foo) <- gene_metadata$ens_id
sce <- sce[rownames(sce) %in% names(foo),]
rownames(sce) <- foo[rownames(sce)]

################
## Parse data ##
################

# names(opts$celltype.colors) <- names(opts$celltype.colors) %>% stringr::str_replace_all("_"," ")
# opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% unique(to.plot$aggregated_celltype)]

###################################
## Plot dimensionality reduction ##
###################################

to.plot <- trajectory.dt %>% merge(sample_metadata.dt, by="cell")

ggplot(to.plot, aes(x=V1, y=V2, fill=celltype)) +
  geom_point(size=1, shape=21, stroke=0.2) +
  scale_fill_manual(values=opts$celltype.colors) +
  ggrepel::geom_text_repel(aes_string(label="celltype"), data=to.plot[,.(V1=median(V1), V2=median(V2)), by="celltype"]) +
  theme_classic() +
  labs(x="", y="") +
  ggplot_theme_NoAxes() +
  theme(
    legend.position="none"
  )


# pdf(paste0(io$outdir,"/rna_umap.pdf"), width=5, height=3, useDingbats = F)
# print(p)
# dev.off()

##########################
## Plot gene expression ##
##########################

genes.to.plot <- c("Dnmt1","Dnmt3a","Dnmt3b","Tet1","Tet2","Tet3")

for (i in genes.to.plot) {
  
  to.plot <- data.table(
    cell = colnames(sce),
    expr = logcounts(sce)[i,]
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
  
  # pdf(sprintf("%s/%s_barplot_pseudobulk.pdf",io$outdir,gene), width=10, height=9)
  print(p)
  # dev.off()
    
}
