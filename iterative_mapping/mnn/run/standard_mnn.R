# Load libraries
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
p$add_argument('--stages',        type="character",   nargs='+',  help='Atlas stage(s)')
p$add_argument('--test_samples',  type="character",   nargs='+',  help='Query batch(es)')
p$add_argument('--test',          action = "store_true",  help = 'Testing mode')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$stages <- c(
#   "E6.5",
#   # "E6.75",
#   # "E7.0",
#   # "E7.25",
#   "E7.5"
#   # "E7.75"
#   # "E8.0",
#   # "E8.25",
#   # "E8.5"
#   # "mixed_gastrulation"
# )
# args$test_samples <- c(
#   # E7.5
#   # "2",
#   # "3",
#   "4",
#   "6"
#   # "19",
#   # "20"
# )
# args$test <- FALSE
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

if (isTRUE(args$test)) print("Test mode activated...")

###############
## Load data ##
###############

source(io$load_data_script)

#############
## Run MNN ##
#############

# Joint normalisation
sce.all <- joint.normalisation(sce.query, sce.atlas, cosineNorm = TRUE)

# Select HVGs
genes <- getHVGs(sce.atlas, block=as.factor(sce.atlas$sample), p.value = 0.10)

# MNN
mapping.dt <- mnn.fn(sce.all, sce.query, sce.atlas, genes = genes, npcs = 50, k = 25)

# table(mapping.dt$celltype_mapped)

##########
## Save ##
##########

fwrite(mapping.dt, sprintf("%s/%s_standard_mnn.txt.gz",io$outdir,paste(args$test_samples,collapse="_")), sep="\t")

