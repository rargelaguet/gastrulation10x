here::i_am("plot_individual_genes/plot_individual_genes.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$outdir <- file.path(io$basedir,"results/individual_genes")

# Options

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

##########
## Plot ##
##########

genes.to.plot <- fread(io$marker_genes)$gene %>% unique %>% .[!grepl("Rik$",.)] %>% head(n=10)
genes.to.plot <- c("Dnmt1","Dnmt3a","Dnmt3b","Tet1","Tet2","Tet3")

genes.to.plot <- genes.to.plot[genes.to.plot%in%rownames(sce)]

for (i in 1:length(genes.to.plot)) {
  
  gene <- genes.to.plot[i]
  
  print(sprintf("%s/%s: %s",i,length(genes.to.plot),gene))
  
  to.plot <- data.table(
    cell = colnames(sce),
    celltype = sce$celltype,
    expr = logcounts(sce)[gene,]
  )
  
  p <- ggplot(to.plot, aes(x=celltype, y=expr, fill=celltype)) +
    # geom_point(shape=21, size=2, data=to.plot[,.(expr=mean(expr)),by="celltype"]) +
    geom_violin(scale="width", alpha=0.75) +
    geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.75) +
    # stat_summary(fun.data = give.n, geom = "text", size=3) +
    scale_fill_manual(values=opts$celltype.colors) +
    labs(x="",y=sprintf("%s expression",gene)) +
    guides(x = guide_axis(angle = 90)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(colour="black",size=rel(0.75)),
      axis.text.y = element_text(colour="black",size=rel(0.9)),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )
  
  pdf(sprintf("%s/%s_boxplots_single_cells.pdf",io$outdir,gene), width=8, height=5)
  print(p)
  dev.off()
    
}


##########
## TEST ##
##########

genes.to.plot <- c("Dppa3")

for (i in genes.to.plot) {
  
  to.plot <- data.table(
    cell = colnames(sce),
    celltype = sce$celltype,
    stage = sce$stage,
    expr = logcounts(sce)[i,]
  ) %>% .[celltype=="PGC"]
  
  p <- ggplot(to.plot, aes(x=stage, y=expr, fill=stage)) +
    # geom_point(shape=21, size=2, data=to.plot[,.(expr=mean(expr)),by="celltype"]) +
    geom_violin(scale="width", alpha=0.75) +
    geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.75) +
    # stat_summary(fun.data = give.n, geom = "text", size=3) +
    scale_fill_manual(values=opts$stage.colors) +
    labs(x="",y=sprintf("%s expression",i)) +
    guides(x = guide_axis(angle = 90)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(colour="black",size=rel(0.75)),
      axis.text.y = element_text(colour="black",size=rel(0.9)),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )
  print(p)
}

