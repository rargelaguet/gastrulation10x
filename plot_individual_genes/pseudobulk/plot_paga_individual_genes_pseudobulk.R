source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$outdir <- paste0(io$basedir,"/results/individual_genes/pseudobulk/paga"); dir.create(io$outdir, showWarnings = F)

# Options
opts$celltypes <- c(
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
	"ExE_endoderm"
	# "Parietal_endoderm",
	# "ExE_ectoderm"
)

opts$rename <- c(
  "Haematoendothelial_progenitors" = "Haemato._progenitors"
)

###############
## Load data ##
###############

# Load SingleCellExperiment object
rna.sce <- readRDS(io$sce_pseudobulk_celltype)[,opts$celltypes]

# Load gene metadata
gene_metadata <- fread(io$gene_metadata) %>%
  .[symbol!="" & ens_id%in%rownames(rna.sce)]

# Rename to gene names
rna.sce <- rna.sce[rownames(rna.sce)%in%gene_metadata$ens_id,]
foo <- gene_metadata$symbol; names(foo) <- gene_metadata$ens_id
new.names <- foo[rownames(rna.sce)]
stopifnot(sum(is.na(new.names))==0)
stopifnot(sum(duplicated(new.names))==0)
rownames(rna.sce) <- new.names

###############
## Load PAGA ##
###############

source(here::here("load_paga_graph.R"))

# Plot graph structure
p <- ggnet2(
  net = net.paga,
  mode = c("x", "y"),
  node.size = 0,
  edge.size = 0.15,
  edge.color = "grey",
  label = FALSE,
  label.size = 2.3
)


###################
## Plot per gene ##
###################

# Define color scale
rna.col.seq <- round(seq(0,1,0.1), 2)
rna.colors <- colorRampPalette(c("gray92", "darkgreen"))(length(rna.col.seq))

# Define genes to plot
# genes.to.plot <- rownames(rna.sce)[grep("Gata",rownames(rna.sce))]
# genes.to.plot <- c("Foxa2","Tfap2a","Mesp1")
# genes.to.plot <- c("Pim2", "Pou5f1", "Slc7a3", "Utf1", "Dppa5a")
genes.to.plot <- c("Dnmt1","Dnmt3a","Dnmt3b","Tet1","Tet2","Tet3")
stopifnot(genes.to.plot%in%rownames(rna.sce))

for (i in 1:length(genes.to.plot)) {
  gene <- genes.to.plot[i]
  print(sprintf("%s/%s: %s",i,length(genes.to.plot),gene))
    
  expr.values <- logcounts(rna.sce[gene,])[1,] %>% minmax.normalisation()
  expr.values[-1] <- expr.values[-1] - 0.05; expr.values[expr.values<0] <- 0
  
  expr.colors <- round(expr.values,1) %>% map(~ rna.colors[which(rna.col.seq == .)]) %>% unlist
  
  p.rna <- p + geom_text(label = "\u25D0", aes(x=x, y=y), color=expr.colors, size=13, family = "Arial Unicode MS",
                      data = p$data[,c("x","y")] %>% dplyr::mutate(expr=expr.colors)) +
    scale_colour_manual(values=expr.colors) + 
    labs(title=gene) +
    theme(
      plot.title = element_text(hjust = 0.5, size=rel(1.50))
    )
  
  
  png(file.path(io$outdir,sprintf("%s_expr_paga.png",gene)), width = 250, height = 280)
  print(p.rna)
  dev.off()
}


#############################
## Plot all genes together ##
#############################
