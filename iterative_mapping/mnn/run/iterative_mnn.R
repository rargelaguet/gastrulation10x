suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--stages',    type="character",   nargs='+',  help='Atlas stage(s)')
p$add_argument('--test_samples',    type="character",   nargs='+',  help='Query batch(es)')
p$add_argument('--test',            action = "store_true",  help = 'Testing mode')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
args$stages <- c(
  "E6.5",
  # "E6.75",
  # "E7.0",
  # "E7.25",
  "E7.5"
  # "E7.75"
  # "E8.0",
  # "E8.25",
  # "E8.5"
  # "mixed_gastrulation"
)
args$test_samples <- c(
  # E7.5
  # "2",
  # "3",
  "4",
  "6"
  # "19",
  # "20"
)
args$test <- FALSE
## END TEST ##

#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation10x/settings.R")
  source("/Users/ricard/gastrulation10x/iterative_mapping/mnn/run/utils.R")
  io$load_data_script <- "/Users/ricard/gastrulation10x/iterative_mapping/mnn/run/load_data.R"
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/gastrulation10x/settings.R")
  source("/homes/ricard/gastrulation10x/iterative_mapping/mnn/run/utils.R")
  io$load_data_script <- "/homes/ricard/gastrulation10x/iterative_mapping/mnn/run/load_data.R"
}  
io$marker_genes <- paste0(io$basedir,"/results/marker_genes/all_stages/marker_genes.txt.gz")
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

# opts$celltypes <- which(table(metadata_atlas$celltype)>25) %>% names
opts$celltypes <- unique(metadata_atlas$celltype)

dist <- fread(io$dist) %>%
  as.data.frame %>% tibble::column_to_rownames("V1") %>% as.matrix %>% as.dist

#######################
## Recursive mapping ##
#######################

sce.query$celltype_mapped <- paste(opts$celltypes,collapse="%")

while (any(grepl("%",sce.query$celltype_mapped))) {
  print(table(sce.query$celltype_mapped))
  mapping.dt <- recursive.fn(sce.query, sce.atlas, dist, cosineNorm = TRUE)
  ids <- match(mapping.dt$cell,colnames(sce.query))
  sce.query$celltype_mapped[ids] <- mapping.dt$celltype_mapped
  sce.query$celltype_score[ids] <- mapping.dt$celltype_score
  sce.query$stage_mapped[ids] <- mapping.dt$stage_mapped
  sce.query$stage_score[ids] <- mapping.dt$stage_score
}

##########
## Save ##
##########

mapping.dt <- data.table(
  cell = colnames(sce.query), 
  celltype_mapped = sce.query$celltype_mapped,
  celltype_score = sce.query$celltype_score,
  stage_mapped = sce.query$stage_mapped,
  stage_score = sce.query$stage_score
)
# foo <- mapping.dt %>% merge(meta_query[,c("cell","celltype.mapped","celltype.score")] %>% setnames(c("cell","celltype_old","score_old")))

fwrite(mapping.dt, sprintf("%s/%s_iterative_mnn.txt.gz",io$outdir,paste(args$test_samples,collapse="_")), sep="\t")
