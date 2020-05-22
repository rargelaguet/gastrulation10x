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
  source("/Users/ricard/gastrulation10x/iterative_mapping/cellassign/utils.R")
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/gastrulation10x/settings.R")
}  
io$outfile <- "/Users/ricard/data/gastrulation10x/results/iterative_mapping/cellassign/test.txt"
dir.create(io$outdir, showWarnings = F)

# Cell types to use
opts$celltypes <- c(
  "Erythroid1", "Erythroid2",
  "Visceral endoderm", "ExE endoderm",
  "Epiblast", "Primitive Streak",
  "ExE ectoderm",
  "Notochord"
)

opts$number.diff.genes <- 50
opts$min.logFC <- 1
opts$threshold.fdr <- 0.01
opts$min_detection_rate_per_group <- 0.25

# Test mode
opts$test_mode <- TRUE

# Maximum of M marker genes per cell type (sorted according to marker score)
# opts$max.genes <- 50

# Update metadata
sample_metadata <- sample_metadata %>%
  .[,celltype:=stringr::str_replace_all(celltype,"_", " ")] %>%
  .[celltype%in%opts$celltypes]

# Subset cells
# if (isTRUE(opts$test_mode)) {
sample_metadata_query <- sample_metadata %>% 
  split(.$celltype) %>% map(~ head(.,n=500)) %>% rbindlist
sample_metadata_atlas <- sample_metadata %>% 
  split(.$celltype) %>% map(~ tail(.,n=1000)) %>% rbindlist
# } else {
#   sample_metadata_query <- sample_metadata
# }

# table(sample_metadata$celltype)

# Load RNA expression data
sce <- readRDS(io$rna.sce)

################
## Parse data ##
################

# remove non-expressed genes
sce <- sce[rowMeans(logcounts(sce))>0,]

# Define query SingleCellExperiment
sce.query <- sce[,sample_metadata_query$cell]
sce.query$celltype <- sample_metadata_query$celltype
dim(sce.query)

# Define atlas SingleCellExperiment
sce.atlas <- sce[,sample_metadata_atlas$cell]
sce.atlas$celltype <- sample_metadata_atlas$celltype
dim(sce.atlas)

########################################################
## Define distance matrix for hierarchical clustering ##
########################################################

io$dist <- "/Users/ricard/data/gastrulation10x/results/phylogenetic_tree/PAGA_distances.csv.gz"
dist <- fread(io$dist) %>%
  as.data.frame %>% tibble::column_to_rownames("V1") %>% as.matrix %>%
  .[opts$celltypes,opts$celltypes] %>% as.dist

plot(hclust(dist))

#######################
## Recursive mapping ##
#######################

sce.query$pred.celltype <- paste(opts$celltypes,collapse="%")

while (any(grepl("%",sce.query$pred.celltype))) {
  print(table(sce.query$pred.celltype))
  dt.pred <- recursive.fn(sce.query, sce.atlas, dist)
  ids <- match(dt.pred$cell,colnames(sce.query))
  sce.query$pred.celltype[ids] <- dt.pred$celltype.pred
}

######################
## Fetch statistics ##
######################

foo <- data.table(cell=colnames(sce.query), celltype.pred = sce.query$pred.celltype)

foobar <- merge(sample_metadata_query[,c("cell","celltype")], foo, by="cell")
mean(foobar$celltype==foobar$celltype.pred)

to.plot <- foobar[,.N, by=c("celltype","celltype.pred")] %>%
  .[,total:=sum(N),by="celltype"] %>% .[,fraction:=N/total] %>% .[,total:=NULL]
#dcast(celltype~celltype.pred, value.var="N")# %>%
# matrix.please

to.plot[,fraction:=round(fraction,2)]

ggplot(to.plot, aes(x=celltype, y=celltype.pred, fill=fraction)) +
  geom_tile(color="black") +
  geom_text(aes(label=fraction)) +
  theme_classic() +
  scale_fill_gradient(low = "gray80", high = "red") +
  labs(x="True", y="Predicted") +
  theme(
    legend.position = "none",
    axis.text.y = element_text(color="black",size=rel(0.8)),
    axis.text.x = element_text(colour="black",size=rel(0.8), angle=90, hjust=1, vjust=0.5)
  )
# pheatmap::pheatmap(to.plot, cluster_rows=F, cluster_cols=F)

# Save cell type predictions
# fwrite(foobar, io$outfile, sep="\t")
