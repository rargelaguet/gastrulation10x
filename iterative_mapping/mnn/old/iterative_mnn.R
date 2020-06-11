suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(argparse))

#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation10x/settings.R")
  source("/Users/ricard/gastrulation10x/iterative_mapping/mnn/utils.R")
  io$load_data_script <- "/Users/ricard/gastrulation10x/iterative_mapping/mnn/run/load_data.R"
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/gastrulation10x/settings.R")
  source("/homes/ricard/gastrulation10x/iterative_mapping/mnn/utils.R")
  io$load_data_script <- "/homes/ricard/gastrulation10x/iterative_mapping/mnn/run/load_data.R"
}  
io$outdir <- paste0(io$basedir,"/results/iterative_mapping/mnn")
io$dist <- paste0(io$basedir,"/results/phylogenetic_tree/PAGA_distances.csv.gz")
dir.create(io$outdir, showWarnings = F)

###############
## Load data ##
###############

source(io$load_data_script)

########################################################
## Define distance matrix for hierarchical clustering ##
########################################################

dist <- fread(io$dist) %>%
  as.data.frame %>% tibble::column_to_rownames("V1") %>% as.matrix %>%
  .[opts$celltypes,opts$celltypes] %>% as.dist

#######################
## Recursive mapping ##
#######################

sce_query$celltype_mapped <- paste(opts$celltypes,collapse="%")

while (any(grepl("%",sce_query$celltype_mapped))) {
  print(table(sce_query$celltype_mapped))
  mapping_dt <- recursive.fn(sce_query, sce_atlas, dist)
  ids <- match(mapping_dt$cell,colnames(sce_query))
  sce_query$celltype_mapped[ids] <- mapping_dt$celltype_mapped
  sce_query$mapping_score[ids] <- mapping_dt$mapping_score
}

##########
## Save ##
##########

mapping_dt <- data.table(
  cell = colnames(sce_query), 
  celltype_mapped = sce_query$celltype_mapped,
  celltype_score = sce_query$celltype_score,
  stage_mapped = sce_query$stage_mapped,
  stage_score = sce_query$stage_score
)
# foo <- mapping_dt %>% merge(meta_query[,c("cell","celltype.mapped","celltype.score")] %>% setnames(c("cell","celltype_old","score_old")))

fwrite(mapping_dt, sprintf("%s/%s_iterative_mnn.txt.gz",io$outdir,paste(args$query_plates,collapse="_")), sep="\t")
