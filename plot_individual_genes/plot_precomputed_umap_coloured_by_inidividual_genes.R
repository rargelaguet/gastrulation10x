here::i_am("plot_individual_genes/plot_precomputed_umap_coloured_by_inidividual_genes.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

io$outdir <- paste0(io$basedir,"/results/individual_genes/umap"); dir.create(io$outdir, showWarnings = F)

opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)

opts$stages <- c(
  "E6.5",
  "E6.75",
  "E7.0",
  "E7.25",
  "E7.5",
  "E7.75",
  "E8.0",
  "E8.25",
  "E8.5",
  "mixed_gastrulation"
)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[stripped==F & doublet==F & celltype%in%opts$celltypes & stage%in%opts$stages] %>%
  .[,celltype:=factor(celltype, levels=opts$celltypes)]

table(sample_metadata$celltype)
table(sample_metadata$stage)


###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- load_SingleCellExperiment(io$sce, cells=sample_metadata$cell, normalise = TRUE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

########################
## Load gene metadata ##
########################

gene_metadata <- fread(io$gene_metadata) %>%
  .[ens_id%in%rownames(sce) & symbol!=""]

##################
## Rename genes ##
##################

foo <- gene_metadata$symbol
names(foo) <- gene_metadata$ens_id
sce <- sce[rownames(sce) %in% names(foo),]
rownames(sce) <- foo[rownames(sce)]

#################################
## Subset SingleCellExperiment ##
#################################

genes.to.plot <- c("Dnmt1","Dnmt3a","Dnmt3b","Tet1","Tet2","Tet3")

# Subset genes
sce_filt <- sce[genes.to.plot,]

# Subset cells
sce_filt <- sce_filt[,sample.int(n=ncol(sce), size=ncol(sce)/5)]

#############################
## Plot one gene at a time ##
#############################

# Extract UMAP coordinates from the metadata
umap.dt <- sample_metadata[,c("umapX","umapY")]

for (i in 1:length(genes.to.plot)) {
  gene <- genes.to.plot[i]
  print(sprintf("%s/%s: %s",i,length(genes.to.plot),gene))
  
  # Create data.table to plot
  to.plot <- data.table(
    cell = colnames(sce),
    expr = logcounts(sce[i,])[1,]
  ) %>% cbind(umap.dt)

  # Plot
  p <- ggplot(to.plot, aes(x=umapX, y=umapY, color=expr)) +
    scale_color_gradient(low = "gray80", high = "red") +
    geom_point(size=0.25) +
    theme_classic() +
    theme(
      legend.title = element_blank(),
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )
    
  # pdf(sprintf("%s/%s.pdf",io$outdir,i), width=8, height=3.5, useDingbats = F)
  # ggsave("ggtest.png", width = 3.25, height = 3.25, dpi = 1200)
  jpeg(sprintf("%s/%s_umap.jpeg",io$outdir,gene), width = 600, height = 600)
  print(p)
  dev.off()
}

#############
## Denoise ##
#############

# NOTE THAT THIS NEEDS A LOT OF MEMORY (AROUND 40GB for 25k CELLS)

pca.mtx <- readRDS(io$pca)[["all"]][colnames(sce_filt),]
rna_denoised.mtx <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(sce_filt)), D=pdist(pca.mtx), k=25)
colnames(rna_denoised.mtx) <- colnames(sce_filt)
assay(sce_filt,"logcounts_denoised") <- rna_denoised.mtx

rm(rna_denoised.mtx); gc(reset=T)

#####################################
## Plot all genes at the same time ##
#####################################

genes.to.plot <- c("Dnmt1","Dnmt3a","Dnmt3b","Tet1","Tet2","Tet3")

to.plot <- as.matrix(assay(sce_filt,"logcounts_denoised")) %>% as.data.table(keep.rownames = T) %>% setnames("rn","gene") %>%
  melt(id.vars="gene", variable.name="cell", value.name="expr") %>%
  merge(sample_metadata[,c("cell","celltype","umapX","umapY")]) %>%
  .[,celltype:=factor(celltype,levels=opts$celltypes)]

to.plot[,expr:=minmax.normalisation(expr),by="gene"]

# Plot
p <- ggplot(to.plot, aes(x=umapX, y=umapY, color=expr)) +
  scale_color_gradient(low = "gray80", high = "red") +
  ggrastr::geom_point_rast(size=0.20) +
  theme_classic() +
  facet_wrap(~gene, ncol=2) +
  ggplot_theme_NoAxes() +
  theme(
    # panel.spacing = unit(2, "lines"),
    strip.text = element_text(colour="black",size=rel(1.25)),
    strip.background = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )

pdf(file.path(io$outdir,"umap_coloured_by_selected_gene_expr.pdf"), width=5.5, height=10)
# jpeg(sprintf("%s/%s_umap.jpeg",io$outdir,gene), width = 600, height = 600)
print(p)
dev.off()
