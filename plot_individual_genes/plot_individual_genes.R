here::i_am("plot_individual_genes/plot_individual_genes.R")

source(here::here("settings.R"))

#####################
## Define settings ##
#####################

## I/O ##

io$outdir <- paste0(io$basedir,"/results/individual_genes")

## Define options ##

# Define stages to plot
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

# Define cell types to plot
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

# Update sample metadata
sample_metadata <- sample_metadata %>% 
  .[stage%in%opts$stages & celltype%in%opts$celltypes] %>%
  .[,celltype:=factor(celltype, levels=opts$celltypes)]

table(sample_metadata$stage)
table(sample_metadata$celltype)

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- readRDS(io$sce)[,sample_metadata$cell]

# Load gene metadata
gene_metadata <- fread(io$gene_metadata) %>%
  .[ens_id%in%rownames(sce)]

# Rename genes
foo <- gene_metadata$symbol
names(foo) <- gene_metadata$ens_id
sce <- sce[rownames(sce) %in% names(foo),]
rownames(sce) <- foo[rownames(sce)]

##########
## Plot ##
##########

# genes.to.plot <- c("Dnmt3l","Apoe")
genes.to.plot <- rownames(sce)

for (i in 1:length(genes.to.plot)) {
  gene <- genes.to.plot[i]
  print(sprintf("%s/%s: %s",i,length(genes.to.plot),gene))
  
  # Create data.table to plot
  to.plot <- data.table(
    cell = colnames(sce),
    celltype = sample_metadata$celltype,
    expr = logcounts(sce[gene,])[1,]
  )

  # Plot
  p <- ggplot(to.plot, aes(x=celltype, y=expr, fill=celltype)) +
    # geom_jitter(size=0.8) +
    geom_violin(scale = "width") +
    geom_boxplot(width=0.5, outlier.shape=NA) +
    scale_fill_manual(values=opts$celltype.colors.1) +
    theme_classic() +
    labs(x="",y="RNA expression") +
    theme(
      axis.text.x = element_text(colour="black",size=rel(1.4), angle=50, hjust=1),
      axis.text.y = element_text(colour="black",size=rel(1.0)),
      legend.position="none"
    )
    
  # pdf(sprintf("%s/%s.pdf",io$outdir,i), width=8, height=3.5, useDingbats = F)
  # ggsave("ggtest.png", width = 3.25, height = 3.25, dpi = 1200)
  jpeg(sprintf("%s/%s.jpeg",io$outdir,gene), width = 900, height = 400)
  print(p)
  dev.off()
}
