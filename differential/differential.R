###########################################################
## Script to do differential expression between lineages ##
###########################################################

suppressMessages(library(SingleCellExperiment))
suppressMessages(library(scater))
suppressMessages(library(edgeR))
suppressMessages(library(argparse))

## Initialize argument parser ##
p <- ArgumentParser(description='')
p$add_argument('--groupA',    type="character",    help='group A')
p$add_argument('--groupB',    type="character",    help='group B')
p$add_argument('--test',      type="character",    help='Statistical test')
p$add_argument('--test_mode', action="store_true", help='Test mode? subset number of cells')
p$add_argument('--outfile',   type="character",    help='Output file')
args <- p$parse_args(commandArgs(TRUE))

# Sanity checks
stopifnot(args$test%in%c("edgeR","t-test","wilcoxon"))

## START TEST
# args$groupA <- c("ExE endoderm")
# args$groupB <- c("Epiblast")
# args$outfile <- c("/Users/ricard/data/gastrulation10x/results/differential/foo.tsv.gz")
# args$outfile <- c("/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/results/differential/foo.tsv.gz")
# args$test_mode <- TRUE
## END TEST


#########
## I/O ##
#########

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation10x/settings.R")
  source("/Users/ricard/gastrulation10x/differential/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation10x/settings.R")
  source("/homes/ricard/gastrulation10x/differential/utils.R")
} else {
  stop("Computer not recognised")
}

#############
## Options ##
#############

# Define groups
args$groupA <- args$groupA %>% stringr::str_replace_all(.,"-"," ")
args$groupB <- args$groupB %>% stringr::str_replace_all(.,"-"," ")
opts$groups <- c(args$groupA,args$groupB)

# Define FDR threshold
opts$threshold_fdr <- 0.01

# Define minimum logFC for significance
opts$min.logFC <- 1

# For a given gene, the minimum fraction of cells that must express it in at least one group
opts$min_detection_rate_per_group <- 0.50


###############
## Load data ##
###############

# Update cell metadata
sample_metadata <- sample_metadata %>%
  .[celltype2%in%opts$groups] %>%
  setnames("celltype2","group") %>%
  .[,c("cell","group")]

if (isTRUE(args$test_mode)) {
  print("Testing mode activated")
  sample_metadata <- sample_metadata %>% split(.,.$group) %>% map(~ head(.,n=100)) %>% rbindlist
}
table(sample_metadata$group)

# Load SingleCellExperiment object
sce <- readRDS(io$rna.sce)[,sample_metadata$cell]
sce$group <- sample_metadata$group

# Load gene metadata
gene_metadata <- fread(io$gene_metadata) %>%
  .[ens_id%in%rownames(sce)] %>%
  .[,c("symbol","ens_id")] %>% 
  setnames("symbol","gene")

# Filter genes
sce <- sce[rownames(sce)%in%gene_metadata$ens_id,]

################################################
## Differential expression testing with edgeR ##
################################################

out <- doDiffExpr(sce, opts$groups, args$test, opts$min_detection_rate_per_group) %>%
  merge(gene_metadata, all.y=T, by="ens_id") %>%
 .[, sig := (padj_fdr<=opts$threshold_fdr & abs(logFC)>opts$min.logFC)] %>%
  setorderv("padj_fdr", na.last=T)

##################
## Save results ##
##################

args$outfile <- args$outfile %>% stringr::str_replace_all(.,"-"," ")
fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
