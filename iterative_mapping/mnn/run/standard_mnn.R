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
#   "E7.5",
#   # "E7.75"
#   # "E8.0",
#   # "E8.25",
#   "E8.5"
#   # "mixed_gastrulation"
# )
# 
# args$test_samples <- c(
#   # E7.5
#   # "2",
#   # "3",
#   "4",
#   "6"
#   # "19",
#   # "20"
# )
# 
# args$test <- FALSE
## END TEST ##

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
io$outdir <- paste0(io$basedir,"/results/iterative_mapping/mnn")
io$dist <- paste0(io$basedir,"/results/phylogenetic_tree/PAGA_distances.csv.gz")

io$outdir <- paste0(io$basedir,"/rna/results/iterative_mapping")

if (isTRUE(args$test)) print("Test mode activated...")

###############
## Load data ##
###############

io$basedir <- "/hps/nobackup2/research/stegle/users/ricard/scnmt_gastrulation_TetChimera/TO_MERGE_scnmt_gastrulation_TetKO"
io$metadata <- paste0(io$basedir,"/sample_metadata.txt.gz")
io$sce <- paste0(io$basedir,"/rna/SingleCellExperiment.rds")

source(io$script_load_data)

#############
## Run MNN ##
#############

mapping_dt <- mnn.fn(sce_query, sce_atlas, npcs = 50, k = 25, cosineNorm = TRUE)

# table(mapping_dt$celltype_mapped)

# foo <- fread("/Users/ricard/data/scnmt_gastrulation_TetChimera/rna/results/iterative_mapping/tet_chimera_march20_plate1_standard_mnn.txt.gz") %>%
#   merge(mapping_dt, by="cell")

# foo <- fread("/Users/ricard/data/scnmt_gastrulation_TetChimera/TO_MERGE_scnmt_gastrulation_TetKO/rna/mapping_10x/sample_metadata_mapping_mnn.txt") %>%
#   merge(mapping_dt %>% copy %>% setnames("cell","id_rna"), by="id_rna") %>%
#   .[,c("id_rna","lineage10x","celltype_mapped")]

##########
## Save ##
##########

fwrite(mapping_dt, sprintf("%s/%s_standard_mnn.txt.gz",io$outdir,paste(args$query_plates,collapse="_")), sep="\t")

