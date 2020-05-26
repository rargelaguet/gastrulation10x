library(ggpubr)

##############
## Settings ##
##############

source("/Users/ricard/gastrulation10x/settings.R")
io$outdir <- paste0(io$basedir,"/results/celltype_affinity"); dir.create(io$outdir, showWarnings = F)

opts$stages <- c(
  "E6.5",
  "E6.75",
  "E7.0",
  "E7.25",
  "E7.5"
  # "E7.75"
  # "E8.0",
  # "E8.25",
  # "E8.5"
  # "mixed_gastrulation"
)

###############
## Load data ##
###############

# Load metadata
sample_metadata <- sample_metadata %>% 
  .[stage%in%opts$stages] %>%
  .[,c("cell","celltype")]

foo <- table(sample_metadata$celltype)<100
if (any(foo)) {
  warning("There are cell types with very small amount of cells, removing them:")
  warning(paste(names(which(foo)),collapse=",  "))
  sample_metadata <- sample_metadata[!celltype%in%names(which(foo))]
}

# Load gene markers
marker_genes.dt <- fread(io$marker_genes) %>%
  .[celltype%in%unique(sample_metadata$celltype)]

# Load SingleCellExperiment
sce <- readRDS(io$rna.sce)[,sample_metadata$cell]

# Add cell types to the SingleCellExperiment object
sce$celltype.pred <- sample_metadata$celltype

##################
## Computations ##
##################

# Calculate atlas-based cell type "score" for each predicted cell type in the query
dt <- unique(sce$celltype.pred) %>% map(function(i) {
# dt <- unique(sce$celltype.pred) %>% head(n=5) %>% map(function(i) {
  logcounts(sce[,sce$celltype.pred==i]) %>% as.matrix %>% 
    rowMeans %>% as.data.table(keep.rownames = T) %>% 
    setnames(c("ens_id","expr")) %>%
    merge(marker_genes.dt[,c("celltype","ens_id")], by="ens_id", allow.cartesian=TRUE) %>%
    .[,.(score=mean(expr)),by="celltype"] %>%
    .[,celltype.pred:=i]
}) %>% rbindlist

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

fwrite(dt, paste0(io$outdir,"/celltype_affinity_E6.5_to_E7.5.txt.gz"))
