here::i_am("plot_individual_genes/pseudobulk/plot_individual_genes_pseudobulk.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$outdir <- file.path(io$basedir,"results/individual_genes/pseudobulk")
# io$sce.pseudobulk <- file.path(io$basedir,"results/pseudobulk/SingleCellExperiment_pseudobulk_celltype.rds")
io$sce.pseudobulk <- file.path(io$basedir,"results/pseudobulk/SingleCellExperiment_pseudobulk_sample_celltype.rds")

# Define options
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

opts$min.cells <- 30

###############
## Load data ##
###############

sce <- readRDS(io$sce.pseudobulk)

# Add metadata
sce$sample <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(1)
sce$celltype <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(2)

# Filter by minimum number of cells
metadata(sce)$n_cells <- metadata(sce)$n_cells[metadata(sce)$n_cells>=opts$min.cells]
sce <- sce[,names(metadata(sce)$n_cells)]

# Filter classes and celltypes manually
sce <- sce[,sce$celltype%in%opts$celltypes]
metadata(sce)$n_cells <- metadata(sce)$n_cells[stringr::str_split(names(metadata(sce)$n_cells), pattern = "-") %>% map_chr(2) %in% opts$celltypes]

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

##############
## Barplots ##
##############

# genes.to.plot <- c("Tex19.1","Morc1","Dppa3","Rex1","Dppa5a","Dppa4","Dppa2","Zfp981")
# genes.to.plot <- c("Pou5f1","Epcam","Fst","Lefty2","Cdkn1c","Acta1")
# genes.to.plot <- grep("^Hox",rownames(sce),value=T)
genes.to.plot <- fread(io$marker_genes)[,gene] %>% unique# %>% head(n=10)

genes.to.plot <- c("Dnmt1","Dnmt3a","Dnmt3b","Tet1","Tet2","Tet3")

genes.to.plot <- genes.to.plot[genes.to.plot%in%rownames(sce)]

for (i in 1:length(genes.to.plot)) {
  
  gene <- genes.to.plot[i]
  print(sprintf("%s/%s: %s",i,length(genes.to.plot),gene))

  to.plot <- data.table(
    expr = logcounts(sce)[gene,],
    celltype = sce$celltype,
    sample = sce$sample
  ) %>% 
    .[,celltype:=factor(celltype,levels=opts$celltypes)]
  
  to.plot[expr==0,expr:=0.10]
  to.plot.means <- to.plot[,.(expr=mean(expr)), by=c("celltype")]
  
  p <- ggplot(to.plot.means, aes(x=celltype, y=expr, fill=celltype)) +
    geom_bar(stat="identity", color="black", width=0.65) +
    geom_jitter(size=1.25, alpha=0.75, width=0.15, shape=21, data=to.plot) +
    # scale_fill_brewer(palette="Dark2") +
    scale_fill_manual(values=opts$celltype.colors) +
    theme_classic() +
    labs(x="",y=sprintf("%s expression",gene)) +
    guides(x = guide_axis(angle = 90)) +
    theme(
      plot.title = element_text(hjust = 0.5),
      strip.text = element_text(size=rel(0.85)),
      axis.text.x = element_text(colour="black",size=rel(0.9)),
      axis.text.y = element_text(colour="black",size=rel(0.9)),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(colour="black",size=rel(1.0)),
      legend.position = "none"
      # legend.title = element_blank(),
      # legend.text = element_text(size=rel(0.85))
    )
  
  pdf(sprintf("%s/%s_barplot_pseudobulk.pdf",io$outdir,gene), width=11, height=6)
  print(p)
  dev.off()
}


