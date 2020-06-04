matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}

# Load libraries
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(edgeR))

#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation10x/settings.R")
  source("/Users/ricard/gastrulation10x/iterative_mapping/mnn/utils.R")
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/gastrulation10x/settings.R")
  source("/homes/ricard/gastrulation10x/iterative_mapping/mnn/utils.R")
}  
io$outfile <- paste0(io$basedir,"/results/iterative_mapping/mnn/test.txt")
io$dist <- paste0(io$basedir,"/results/phylogenetic_tree/PAGA_distances.csv.gz")
dir.create(io$outdir, showWarnings = F)

# Cell types to use
opts$celltypes <- c(
  "Erythroid1", "Erythroid2",
  "Visceral endoderm", "ExE endoderm",
  "Epiblast", "Primitive Streak",
  "ExE ectoderm",
  "Notochord"
)
# opts$celltypes <- opts$celltypes %>% stringr::str_replace_all("/", "_") %>% stringr::str_replace_all("_", " ") 

# Test mode
opts$test_mode <- TRUE

# Update metadata
sample_metadata <- sample_metadata %>%
  .[,celltype:=stringr::str_replace_all(celltype,"/", "_")] %>%
  .[,celltype:=stringr::str_replace_all(celltype,"_", " ")] %>%
  .[celltype%in%opts$celltypes]

sample_metadata_atlas <- sample_metadata %>% .[stage%in%c("E7.0","E7.25","E7.75","E8.0")]
sample_metadata_query <- sample_metadata %>% .[stage=="E7.5"]
  
# Subset cells
if (isTRUE(opts$test_mode)) {
  sample_metadata_query <- sample_metadata %>% split(.$celltype) %>% map(~ head(.,n=50)) %>% rbindlist
  sample_metadata_atlas <- sample_metadata %>% split(.$celltype) %>% map(~ tail(.,n=50)) %>% rbindlist
}

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

dist <- fread(io$dist) %>%
  as.data.frame %>% tibble::column_to_rownames("V1") %>% as.matrix %>%
  .[opts$celltypes,opts$celltypes] %>% as.dist

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

pred.dt <- data.table(cell=colnames(sce.query), celltype.pred = sce.query$pred.celltype) %>%
  merge(sample_metadata_query[,c("cell","celltype")])
mean(pred.dt$celltype==pred.dt$celltype.pred)

##########
## Plot ##
##########

to.plot <- pred.dt[,.N, by=c("celltype","celltype.pred")] %>%
  .[,total:=sum(N),by="celltype"] %>% .[,fraction:=N/total] %>% .[,total:=NULL]

p <- ggplot(to.plot, aes(x=celltype, y=celltype.pred, fill=fraction)) +
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

pdf(paste0(io$outdir,"/foo.pdf"), width=8, height=6.5, useDingbats = F)
print(p)
dev.off()

##########
## Save ##
##########

fwrite(pred.dt, io$outfile, sep="\t")
