library(ggpubr)

##############
## Settings ##
##############

source("/Users/ricard/gastrulation10x/settings.R")
io$outdir <- paste0(io$basedir,"/results/celltype_score"); dir.create(io$outdir, showWarnings = F)

opts$scale <- FALSE

###############
## Load data ##
###############

# Load metadata
sample_metadata <- sample_metadata %>% 
  .[,c("cell","celltype")]

# Load gene markers
marker_genes.dt <- fread(io$marker_genes)

# Load SingleCellExperiment
sce <- readRDS(io$rna.sce)[,sample_metadata$cell]

# Add predicted cell types to the SCE object
sce$celltype.pred <- sample_metadata$celltype

################
## Parse data ##
################

if (isTRUE(opts$scale)) {
  stop("TO-DO")
} else {
  expr.matrix <- logcounts(sce)
}

##################
## Computations ##
##################

# Calculate atlas-based cell type "score" for each predicted cell type in the query
dt <- unique(sce$celltype.pred) %>% map(function(i) {
# dt <- unique(sce$celltype.pred) %>% head(n=5) %>% map(function(i) {
  expr.matrix[,sce$celltype.pred==i] %>% as.matrix %>% 
    rowMeans %>% as.data.table(keep.rownames = T) %>% 
    setnames(c("ens_id","expr")) %>%
    merge(marker_genes.dt[,c("celltype","ens_id")], by="ens_id", allow.cartesian=TRUE) %>%
    .[,.(score=mean(expr)),by="celltype"] %>%
    .[,celltype.pred:=i]
}) %>% rbindlist

# TO-DO Calculate cell type "score" for each single cell

##########
## Plot ##
##########

colors <- opts$celltype.colors.1
names(colors) <- names(colors) %>% stringr::str_replace_all("_"," ")

# Plot number of marker genes per cell types
for (i in unique(dt$celltype.pred)) {
  
  to.plot <- dt[celltype.pred==i] %>%
    .[,celltype:=stringr::str_replace_all(celltype,"_"," ")] %>%
    .[,celltype:=factor(celltype,levels=names(colors))]
  
  p <- ggbarplot(to.plot, x="celltype", y="score", fill="celltype") +
    scale_fill_manual(values=colors) +
    labs(x="", y="Cell type affinity") +
    theme(
      axis.text.y = element_text(size=rel(0.75)),
      axis.text.x = element_text(colour="black",size=rel(0.8), angle=90, hjust=1, vjust=0.5),
      axis.ticks.x = element_blank(),
      legend.position = "none"
  )
  
  pdf(sprintf("%s/%s_affinity.pdf",io$outdir,i), width = 9, height = 5)
  print(p)
  dev.off()
}

##########
## Save ##
##########

# fwrite(dt.filt, paste0(io$outdir,"/marker_genes.txt.gz"))
