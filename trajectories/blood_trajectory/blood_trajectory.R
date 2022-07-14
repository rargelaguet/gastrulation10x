here::i_am("trajectories/blood_trajectory/blood_trajectory.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## I/O and options ##
#####################

io$outdir <- paste0(io$basedir,"/results/trajectories/blood_precomputed")

names(opts$celltype.colors) <- names(opts$celltype.colors) %>% stringr::str_replace_all("_", " ")

opts$haem_colours = c(
  "Mes1"= "#c4a6b2",
  "Mes2"= "#ca728c",
  
  "Cardiomyocytes" =  "#B51D8D",  
  
  "BP1" = "#6460c5",
  "BP2" = "#96b8e4",
  "Haem3"= "#02f9ff",
  "BP3" = "#07499f",
  "BP4" = "#036ef9",
  
  "Haem1"= "#bb22a7",
  "Haem2" = "#f695e9",
  "Haem4" = "#4c4a81",
  
  "EC1"= "#006737",
  
  "EC2" = "#5eba46",
  "EC3" = "#818068",
  "EC4"="#d6de22",
  "EC5"="#5f7238",
  "EC6"="#2ab572",
  "EC7"="#000000",
  "EC8"="#a0cb3b",
  
  "Ery1"="#f67a58",
  "Ery2" ="#a26724",
  "Ery3"="#cdaf7f",
  "Ery4"= "#625218",
  
  "My" = "#c62127",
  "Mk"= "#f6931d"
)

#####################
## Load trajectory ##
#####################

trajectory.dt <- fread(file.path(io$basedir,"original/sample_metadata.txt.gz")) %>%
  .[!is.na(haem_gephiY)] %>%
  .[,c("cell","sample", "stage", "celltype", "umapX", "umapY", "haem_gephiX", "haem_gephiY","haem_subclust")] %>%
  setnames(c("haem_gephiX","haem_gephiY"),c("V1","V2")) 

for (i in c("umapX", "umapY", "V1", "V1")) {
  trajectory.dt[[i]] <- round(trajectory.dt[[i]],2)
}

# fwrite(trajectory.dt, file.path(io$outdir,"blood_trajectory.txt.gz"), sep="\t")

#########################
## Load RNA expression ##
#########################

sce <- load_SingleCellExperiment(io$sce, normalise = T, cells=trajectory.dt$cell)

# Add sample metadata as colData
colData(sce) <- trajectory.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

# Rename genes
gene_metadata <- fread(io$gene_metadata) %>% .[ens_id%in%rownames(sce)]
foo <- gene_metadata$symbol
names(foo) <- gene_metadata$ens_id
sce <- sce[rownames(sce) %in% names(foo),]
rownames(sce) <- foo[rownames(sce)]

# Remove lowly expressed genes
# tmp <- rowSums(counts(sce))
# sce <- sce[tmp>=25,]

# Save
# saveRDS(sce,file.path(io$outdir,"SingleCellExperiment.rds"))

##########################################################
## Plot dimensionality reduction coloured per cell type ##
##########################################################

to.plot <- trajectory.dt

p <- ggplot(to.plot, aes(x=V1, y=V2)) +
  ggrastr::geom_point_rast(aes(fill=celltype), size=1, shape=21, stroke=0.1, raster.dpi=150) + 
  scale_fill_manual(values=opts$celltype.colors) +
  theme_classic() +
  ggplot_theme_NoAxes() +
  theme(
    legend.position = "none",
  )

# Label clusters
groups <- unique(to.plot$celltype)
labels.loc <- lapply(groups, function(i) {
  data.table(t(apply(X = to.plot[celltype==i, c("V1","V2")], MARGIN = 2, FUN = median, na.rm = TRUE))) %>%
    .[,celltype:=i] %>% return
}) %>% rbindlist %>% setnames(c("V1","V2","celltype"))
p <- p + geom_text_repel(aes_string(label="celltype"), data=labels.loc)

pdf(paste0(io$outdir,"/blood_trajectory_celltype.pdf"), width=6, height=6)
print(p)
dev.off()

#################################################################
## Plot dimensionality reduction coloured per blood subcluster ##
#################################################################

to.plot <- trajectory.dt

p <- ggplot(to.plot, aes(x=V1, y=V2)) +
  ggrastr::geom_point_rast(aes(colour=haem_subclust), size=0.05) + scale_color_manual(values=opts$haem_colours) +
  # ggrastr::geom_point_rast(aes(colour=celltype), size=0.05) + scale_color_manual(values=opts$celltype.colors) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )

# Label clusters
groups <- unique(to.plot$haem_subclust)
labels.loc <- lapply(groups, function(i) {
    data.table(t(apply(X = to.plot[haem_subclust==i, c("V1","V2")], MARGIN = 2, FUN = median, na.rm = TRUE))) %>%
      .[,haem_subclust:=i] %>% return
  }) %>% rbindlist %>% setnames(c("V1","V2","haem_subclust"))
p <- p + geom_text_repel(aes_string(label="haem_subclust"), data=labels.loc)

pdf(paste0(io$outdir,"/blood_trajectory_haem_subclust.pdf"), width=7, height=6)
print(p)
dev.off()


####################################
## Denoise RNA expression for viz ##
####################################

# Select genes
# genes.to.plot <- c("Dnmt1","Dnmt3a","Dnmt3b","Tet1","Tet2","Tet3")
genes.to.plot <- c("Hbb-y","Hbb-bh1","Hba-x","Hba-a1","Hba-a2","Runx1","Gata1","Klf1")
# genes.to.plot <- grep("Hb",rownames(sce),value=T)
sce_filt <- sce[genes.to.plot]

trajectory.mtx <- trajectory.dt[,c("cell","V1","V2")] %>% matrix.please %>% .[colnames(sce_filt),]
rna_denoised.mtx <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(sce_filt)), D=pdist(trajectory.mtx), k=25)
colnames(rna_denoised.mtx) <- colnames(sce_filt)
assay(sce,"logcounts_denoised") <- rna_denoised.mtx
rm(rna_denoised.mtx); gc(reset=T)

#########################
## Plot RNA expression ##
#########################

# Subset cells
# sce_filt <- sce_filt[,sample.int(n=ncol(sce_filt), size=ncol(sce_filt)/1.5)]

# Plot
for (i in genes.to.plot) {
  
  to.plot <- data.table(
    cell = colnames(sce_filt),
    expr = logcounts(sce_filt)[i,]
    # expr_denoised = assay(sce_filt,"logcounts_denoised")[i,]
  ) %>% 
    melt(id.vars=c("cell"), value.name="expr") %>%
    merge(trajectory.dt, by="cell")
  
  p <- ggplot(to.plot[variable=="expr"], aes(x=V1, y=V2)) +
    ggrastr::geom_point_rast(color="gray90", size=0.75, alpha=0.5, data=to.plot[expr==0]) +
    ggrastr::geom_point_rast(aes(fill=expr), size=1, shape=21, stroke=0.02, data=to.plot[expr>0]) +
    scale_fill_gradient2(low = "gray50", mid="gray90", high = "darkred") +
    # facet_wrap(~variable, nrow=1) +
    # ggrepel::geom_text_repel(aes_string(label="celltype"), data=to.plot[,.(V1=median(V1), V2=median(V2)), by="celltype"]) +
    theme_classic() +
    labs(x="", y="", title=i) +
    ggplot_theme_NoAxes() +
    theme(
      legend.position="none"
    )
  
  pdf(sprintf("%s/blood_trajectory_expr_%s.pdf",io$outdir,i), width=5, height=5)
  print(p)
  dev.off()
  
}
