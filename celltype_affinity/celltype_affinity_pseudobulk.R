
##############
## Settings ##
##############

source("/Users/ricard/gastrulation10x/settings.R")
io$outdir <- paste0(io$basedir,"/results/celltype_affinity")

opts$stages <- c(
  "E6.5",
  "E6.75",
  "E7.0",
  "E7.25",
  "E7.5",
  "E7.75",
  "E8.0",
  "E8.25",
  "E8.5"
  # "mixed_gastrulation"
)

#####################
## Update metadata ##
#####################

sample_metadata <- sample_metadata %>% 
  .[stage%in%opts$stages] %>%
  .[,c("cell","celltype")]

# foo <- table(sample_metadata$celltype)<100
# if (any(foo)) {
#   warning("There are cell types with very small amount of cells, removing them:")
#   warning(paste(names(which(foo)),collapse=",  "))
#   sample_metadata <- sample_metadata[!celltype%in%names(which(foo))]
# }

###############
## Load data ##
###############

# Load gene markers
marker_genes.dt <- fread(io$marker_genes) %>%
  .[celltype%in%unique(sample_metadata$celltype)]
length(unique(marker_genes.dt$celltype))

# Load SingleCellExperiment
# sce <- readRDS(io$rna.sce)[,sample_metadata$cell]

# Load average expression per celltype and gene
pseudobulk_expr.dt <- fread(io$average_expression_per_celltype) %>%
  .[ens_id%in%unique(marker_genes.dt$ens_id) & gene!=""] %>%
  .[,group:=stringr::str_replace_all(group," ", "_")] %>%
  .[group=="Forebrain/Midbrain/Hindbrain",group:="Forebrain_Midbrain_Hindbrain"] %>%
  setnames("group","celltype") %>%
  .[celltype%in%unique(sample_metadata$celltype)]

##################
## Computations ##
##################

# Calculate correlation coefficient between each pair of cell types, across genes 
m <- pseudobulk_expr.dt %>% 
  .[,id:=paste(ens_id,gene,sep="_")] %>%
  dcast(id~celltype, value.var="mean_expr") %>%
  matrix.please

r <- cor(m)
# diag(r) <- NA

##########
## Plot ##
##########

pdf(sprintf("%s/correlation_celltypes_heatmap.pdf",io$outdir), width = 11, height = 8)
pheatmap::pheatmap(r)
dev.off()

