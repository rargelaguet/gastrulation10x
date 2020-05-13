matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}

# Set up reticulate connection
suppressPackageStartupMessages(library(reticulate))
if (grepl("ricard",Sys.info()['nodename'])) {
  use_python("/Users/ricard/anaconda3/envs/gpflow_v2/bin/python", required = TRUE)
} else if(grepl("ebi",Sys.info()['nodename'])){
  use_python("/nfs/research1/stegle/users/ricard/conda-envs/gpflow2/bin/python", required = TRUE)
}  
py_config()

# Load libraries
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tensorflow))
suppressPackageStartupMessages(library(cellassign))
suppressPackageStartupMessages(library(edgeR))

#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation10x/settings.R")
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/gastrulation10x/settings.R")
}  
io$outdir <- "/Users/ricard/data/gastrulation10x/results/cellassign/test"
dir.create(io$outdir, showWarnings = F)

# opts$stages <- c("E6.5")

# Cell types to use
# opts$celltypes <- opts$celltypes.1
opts$celltypes <- c("Erythroid1","Erythroid2","Visceral endoderm", "ExE endoderm")

# Maximum of M marker genes per cell type (sorted according to marker score)
opts$max.genes <- 50

# Update metadata
sample_metadata <- sample_metadata %>%
.[celltype%in%opts$celltypes]

# Subset cells
sample_metadata <- sample_metadata %>% split(.$celltype) %>% 
    map(~ head(.,n=100)) %>% rbindlist

###############
## Load data ##
###############

# Load hierarchical clustering
h <- readRDS("/Users/ricard/data/gastrulation10x/results/phylogenetic_tree/hclust.rds")

# Load RNA expression data
sce <- readRDS(io$rna.sce)


#############################
## Hierarchical clustering ##
#############################

dt <- fread("/Users/ricard/data/gastrulation10x/results/marker_genes/avg_expr_per_celltype_and_gene.txt.gz") %>%
  .[,gene:=NULL] %>% unique

marker_genes <- fread("/Users/ricard/data/gastrulation10x/results/marker_genes/marker_genes.txt.gz") %>%
  .[,head(.SD,n=50),by="celltype"]

m <- dt %>% 
  .[ens_id%in%unique(marker_genes$ens_id)] %>%
  dcast(ens_id~group,value.var="mean_expr") %>% matrix.please %>% t

h <- hclust(dist(m))

plot(h)

#######################
## Recursive mapping ##
#######################

recursive_mapping.fn <- function(sce, h) {
  
  # Hierarchical clustering
  
  # Cut the tree
  cut <- cutree(h,k=2)
  
  # Define groups
  groupA <- names(which(cut==1))
  groupB <- names(which(cut==2))
  
  if (length(groupA)==1 & length(groupB)==1) {
    return()
  }
  stopifnot(length(intersect(groupA,groupB))==0)
  
  sample_metadata.filt <- sample_metadata %>%
    .[celltype%in%c(groupA,groupB)] %>%
    .[,group:=c("groupA","groupB")[as.numeric(celltype%in%groupA)+1]] %>%
    .[,group:=factor(group,levels=c("groupA","groupB"))] %>% 
    setorder(group)
  
  table(sample_metadata.filt$group)
  
  sce.filt <- sce[,sample_metadata.filt$cell]
  sce.filt$group <- sample_metadata.filt$group
  
  fdr_threshold <- 0.01
  
  diff <- differential_expression(sce, groupA, groupB)
  markers.groupA <- diff[padj_fdr<fdr_threshold & logFC>1,ens_id] %>% head(n=25)
  markers.groupB <- diff[padj_fdr<fdr_threshold & logFC<(-1),ens_id] %>% head(n=25)
  
  marker_list <- list(
    "groupA" = markers.groupA,
    "groupB" = markers.groupB
  )
  
  cellassign.fit <- cellassign.fn(sce, marker_list)
  celltypes(cellassign.fit)
    
}



differential_expression <- function(sce.filt) {
  stopifnot(!is.null(sce.filt$group))    
  
  groups <- unique(sce.filt$group)
  
  # Filter genes by detection rate per group
  min_detection_rate_per_group <- 0.25
  cdr_A <- rowMeans(logcounts(sce.filt[,sce.filt$group==groups[1]])>0) >= min_detection_rate_per_group
  cdr_B <- rowMeans(logcounts(sce.filt[,sce.filt$group==groups[2]])>0) >= min_detection_rate_per_group
  sce.filt <- sce.filt[cdr_B | cdr_A,]

  # Convert SCE to DGEList
  sce_edger <- scran::convertTo(sce.filt, type="edgeR")
  
  # Define design matrix
  cdr <- colMeans(logcounts(sce.filt)>0)
  design <- model.matrix(~cdr+sce.filt$group)
  
  # Estimate dispersions
  sce_edger <- estimateDisp(sce_edger,design)
  
  # Fit GLM
  fit <- glmQLFit(sce_edger,design)
  
  # Likelihood ratio test
  lrt <- glmQLFTest(fit)
  
  # Construct output data.frame
  out <- topTags(lrt, n=nrow(lrt))$table %>% 
    as.data.table(keep.rownames=T) %>%
    setnames(c("ens_id","logFC","logCPM","LR","p.value","padj_fdr")) %>%
    .[,log_padj_fdr:= -log10(padj_fdr)] %>%
    .[,c("logCPM","LR"):=NULL] %>%
    setorder(padj_fdr)
  
  
  return(out)
}

cellassign.fn <- function(sce.filt, marker_list) {
  
  # Create binary membership matrix  
  bmat <- marker_list_to_mat(marker_list)
  
  # Subset SingleCellExperiment
  sce.filt <- sce.filt[rownames(bmat),]
  
  # Extract size factors
  s <- sizeFactors(sce.filt)
  
  # Run
  fit <- cellassign(sce.filt, 
    marker_gene_info = bmat, 
    min_delta = 1,
    s = s, 
    # learning_rate = 1e-2, 
    shrinkage = TRUE,
    verbose = FALSE
  )
  return(fit)
}

##################
## Query output ##
##################

# Plot heatmap of cell type probabilities
pheatmap::pheatmap(cellprobs(fit))

# Compare to ground truth
foo <- prop.table(table(sce.filt$group, celltypes(fit)), margin = 1)
pheatmap::pheatmap(foo, cluster_rows = F, cluster_cols = F)
