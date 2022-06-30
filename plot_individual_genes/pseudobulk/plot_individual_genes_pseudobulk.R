here::i_am("plot_individual_genes/pseudobulk/plot_individual_genes_pseudobulk.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$outdir <- file.path(io$basedir,"results/individual_genes/pseudobulk")
io$sce.pseudobulk <- file.path(io$basedir,"results/pseudobulk/SingleCellExperiment_pseudobulk_celltype.rds")

# Define options
opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  # "PGC",
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
  "Visceral_endoderm"
  # "ExE_endoderm",
  # "ExE_ectoderm",
  # "Parietal_endoderm"
)

###############
## Load data ##
###############

sce <- readRDS(io$sce.pseudobulk)[,opts$celltypes]

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

#############################
## Ploe one gene at a time ##
#############################

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
    celltype = colnames(sce),
    expr = logcounts(sce)[gene,]
  )
  
  p <- ggplot(to.plot, aes(x=celltype, y=expr, fill=celltype)) +
    geom_bar(stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors) +
    labs(x="",y=sprintf("%s expression",gene)) +
    guides(x = guide_axis(angle = 90)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(colour="black",size=rel(1)),
      axis.text.y = element_text(colour="black",size=rel(0.9)),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )
  
  pdf(sprintf("%s/%s_barplot_pseudobulk.pdf",io$outdir,gene), width=11, height=6)
  print(p)
  dev.off()
}


#######################################
## Plot multple genes simultaneously ##
#######################################

# genes.to.plot <- c("Dnmt1","Dnmt3a","Dnmt3b","Tet1","Tet2","Tet3")
genes.to.plot <- c("Itga2b","Kdr")

to.plot <- logcounts(sce)[genes.to.plot,] %>% as.data.table(keep.rownames = T) %>% setnames("rn","gene") %>%
   melt(id.vars="gene", variable.name="celltype", value.name="expr") %>%
  .[,celltype:=factor(celltype,levels=opts$celltypes)]

p <- ggplot(to.plot, aes(x=celltype, y=expr, fill=celltype)) +
  geom_bar(stat="identity", color="black") +
  scale_fill_manual(values=opts$celltype.colors[opts$celltypes]) +
  facet_wrap(~gene) +
  labs(x="",y="Gene expression") +
  guides(x = guide_axis(angle = 90)) +
  theme_classic() +
  theme(
    panel.spacing = unit(2, "lines"),
    strip.text = element_text(colour="black",size=rel(1.25)),
    strip.background = element_blank(),
    axis.text.x = element_text(colour="black",size=rel(0.75)),
    axis.text.y = element_text(colour="black",size=rel(0.9)),
    # axis.ticks.x = element_blank(),
    legend.position = "none"
  )


# pdf(sprintf("%s/barplot_pseudobulk_dnmt_tet.pdf",io$outdir), width=14, height=8)
pdf(sprintf("%s/barplot_pseudobulk_dnmt_tet.pdf",io$outdir), width=14, height=8)
print(p)
dev.off()


#################
## Scatterplot ##
#################

gene.x <- "Itga2b"
gene.y <- "Kdr"

to.plot <- data.table(
  celltype = colnames(sce),
  x = logcounts(sce)[gene.x,],
  y = logcounts(sce)[gene.y,]
)
to.plot[,c("x","y"):=list(x/max(x), y/max(y))]

to.plot[,dot_size:=minmax.normalisation(x*y)]
to.plot[,dot_size:=max(x,y),by="celltype"] %>% .[,dot_size:=minmax.normalisation(dot_size)]

p <- ggscatter(to.plot, x="x", y="y", fill="celltype", shape=21, size="dot_size") +
  labs(x=sprintf("%s expression",gene.x),y=sprintf("%s expression",gene.y)) +
  scale_fill_manual(values=opts$celltype.colors[opts$celltypes]) +
  geom_text(aes(label=celltype), data=to.plot[x>=0.65 | y>=0.65], size=3) +
  scale_size_continuous(range = c(2,10)) +
  scale_x_continuous(breaks=c(0,0.5,1)) + scale_y_continuous(breaks=c(0,0.5,1)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(colour="black",size=rel(1)),
    axis.text.y = element_text(colour="black",size=rel(1)),
    legend.position = "none"
  )

pdf(sprintf("%s/scatterplot_%s_vs_%s.pdf",io$outdir,gene.x,gene.y), width=6, height=5)
print(p)
dev.off()
