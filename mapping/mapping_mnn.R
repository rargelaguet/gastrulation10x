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
p$add_argument('--stages',       type="character",   nargs='+',  help='Atlas stage(s)')
p$add_argument('--celltypes',    type="character",   nargs='+',  help='Atlas cell types(s)')
p$add_argument('--test_samples', type="character",   nargs='+',  help='Query batch(es)')
p$add_argument('--nPCs',         type="integer",                 help='Number of PCs')
p$add_argument('--k',            type="integer",                 help='Number of neighbours')
p$add_argument('--test',         action = "store_true",          help = 'Testing mode')
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
# args$celltypes = c(
#   # "Epiblast"
#   # "Primitive_Streak",
#   # "Caudal_epiblast",
#   # "PGC",
#   # "Anterior_Primitive_Streak",
#   # "Notochord",
#   # "Def._endoderm",
#   # "Gut",
#   # "Nascent_mesoderm",
#   # "Mixed_mesoderm",
#   # "Intermediate_mesoderm",
#   # "Caudal_Mesoderm",
#   # "Paraxial_mesoderm",
#   # "Somitic_mesoderm",
#   # "Pharyngeal_mesoderm",
#   # "Cardiomyocytes",
#   # "Allantois",
#   # "ExE_mesoderm",
#   # "Mesenchyme",
#   # "Haematoendothelial_progenitors",
#   # "Endothelium",
#   # "Blood_progenitors_1",
#   # "Blood_progenitors_2",
#   # "Erythroid1",
#   # "Erythroid2",
#   # "Erythroid3",
#   # "NMP",
#   # "Rostral_neurectoderm",
#   # "Caudal_neurectoderm",
#   # "Neural_crest",
#   # "Forebrain_Midbrain_Hindbrain",
#   # "Spinal_cord",
#   # "Surface_ectoderm",
#   # "Visceral_endoderm",
#   # "ExE_endoderm",
#   "ExE_ectoderm"
#   # "Parietal_endoderm"
# )
# args$celltypes <- NULL
# args$test_samples <- c(
#   # E7.5
#   # "2",
#   # "3",
#   "4"
#   # "6"
#   # "19",
#   # "20"
# )
# args$test <- FALSE
# args$nPCs <- 5
# args$k <- 10
## END TEST ##

#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation10x/settings.R")
  source("/Users/ricard/gastrulation10x/mapping/utils.R")
  io$load_data_script <- "/Users/ricard/gastrulation10x/mapping/load_data.R"
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/gastrulation10x/settings.R")
  source("/homes/ricard/gastrulation10x/mapping/utils.R")
  io$load_data_script <- "/homes/ricard/gastrulation10x/mapping/load_data.R"
}  
# io$marker_genes <- paste0(io$basedir,"/results/marker_genes/all_stages/marker_genes.txt.gz")
io$outdir <- paste0(io$basedir,"/results/mapping/stages")

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
sce.atlas <- multiBatchNorm(sce.atlas, batch=as.factor(sce.atlas$sample))
genes <- getHVGs(sce.atlas, p.value = 0.10)
# genes <- getHVGs(sce.atlas, block=as.factor(sce.atlas$sample), p.value = 0.10)
# genes <- getHVGs(sce.all, block=sce.all$block, p.value = 0.10)

# MNN
mapping.dt <- mnn.fn(sce.all, sce.query, sce.atlas, genes = genes, npcs = args$nPCs, k = args$k)

# table(mapping.dt$celltype_mapped)

##########
## Save ##
##########

if (!is.null(args$celltypes)) {
  fwrite(mapping.dt, sprintf("%s/mapping_mnn_%s_%s.txt.gz",io$outdir,paste(args$test_samples,collapse="_"),paste(args$celltypes,collapse="_")), sep="\t")
} else {
  fwrite(mapping.dt, sprintf("%s/mapping_mnn_%s.txt.gz",io$outdir,paste(args$test_samples,collapse="_")), sep="\t")
  
}

