
suppressMessages(library(scater))
suppressMessages(library(edgeR))
suppressMessages(library(argparse))

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation10x/settings.R")
  source("/Users/ricard/gastrulation10x/utils.R")
  source("/Users/ricard/gastrulation10x/differential/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation10x/settings.R")
  source("/homes/ricard/gastrulation10x/utils.R")
  source("/homes/ricard/gastrulation10x/differential/utils.R")
} else {
  stop("Computer not recognised")
}

io$TFs <- paste0(io$basedir,"/results/differential/celltypes/TFs/TFs.txt")

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--groupA',    type="character",    help='group A')
p$add_argument('--groupB',    type="character",    help='group B')
p$add_argument('--test_mode', action="store_true", help='Test mode? subset number of cells')
p$add_argument('--outfile',   type="character",    help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
# args$groupA <- c("Mixed_mesoderm")
# args$groupB <- c("ExE_ectoderm")
# args$outfile <- c("/Users/ricard/data/gastrulation10x/results/differential/TFs/foo.tsv.gz")
# args$outfile <- c("/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/results/differential/foo.tsv.gz")
# args$test_mode <- FALSE
## END TEST

#############
## Options ##
#############

# Sanity checks
stopifnot(args$groupA%in%opts$celltypes)
stopifnot(args$groupB%in%opts$celltypes)

# Define groups
opts$groups <- c(args$groupA,args$groupB)

# Define FDR threshold
opts$threshold_fdr <- 0.01

# Define minimum logFC for significance
opts$min.logFC <- 1

# For a given gene, the minimum fraction of cells that must express it in at least one group
opts$min_detection_rate_per_group <- 0.25

###################
## Load metadata ##
###################

# Update cell metadata
sample_metadata <- sample_metadata %>%
  .[celltype%in%opts$groups] %>%
  setnames("celltype","group") %>%
  .[,c("cell","group")]

# Sort cells so that groupA comes before groupB
sample_metadata[,group:=factor(group,levels=opts$groups)] %>% setorder(group)

if (isTRUE(args$test_mode)) {
  print("Testing mode activated")
  sample_metadata <- sample_metadata %>% split(.,.$group) %>% map(~ head(.,n=250)) %>% rbindlist
}
table(sample_metadata$group)

###############################
## Load gene and TF metadata ##
###############################

gene_metadata <- fread(io$gene_metadata) %>%
  .[,c("symbol","ens_id")] %>% 
  setnames("symbol","gene")

# Load TF information
opts$TFs <- fread(io$TFs)[[1]] %>% stringr::str_to_title()
TF_metadata <- gene_metadata[gene%in%opts$TFs]


######################################
## Load SingleCellExperiment object ##
######################################

# sce <- readRDS(io$rna.sce)[,sample_metadata$cell]
sce <- load_SingleCellExperiment(io$sce, normalise = TRUE, cells = sample_metadata$cell)
sce$group <- sample_metadata$group

# Subset SingleCellExperiment to TFs
sce <- sce[rownames(sce)%in%TF_metadata$ens_id]
TF_metadata <- TF_metadata[ens_id%in%rownames(sce)]

################
## Parse data ##
################

# calculate detection rate per gene
cdr.dt <- data.table(
  ens_id = rownames(sce),
  detection_rate_A = rowMeans(logcounts(sce[,sce$group==opts$groups[1]])>0),
  detection_rate_B = rowMeans(logcounts(sce[,sce$group==opts$groups[2]])>0)
) %>% setnames(c("ens_id",sprintf("detection_rate_%s",opts$groups[1]),sprintf("detection_rate_%s",opts$groups[2])))
# .[,cdr_diff:=abs(out[,(sprintf("detection_rate_%s",opts$groups[1])),with=F][[1]] - out[,(sprintf("detection_rate_%s",opts$groups[2])),with=F][[1]])] %>%

################################################
## Differential expression testing with edgeR ##
################################################

out <- doDiffExpr(sce, groups = opts$groups, test = "edgeR", min_detection_rate_per_group = opts$min_detection_rate_per_group) %>%
  merge(cdr.dt, all.y=T, by="ens_id") %>%
  merge(TF_metadata, all.y=T, by="ens_id") %>%
 .[, sig := (padj_fdr<=opts$threshold_fdr & abs(logFC)>=opts$min.logFC)] %>%
  # setorderv(c("sig","padj_fdr"), na.last=T)
  setorder(-sig, padj_fdr, na.last=T)

##################
## Save results ##
##################

# args$outfile <- args$outfile# %>% stringr::str_replace_all(.,"-"," ")
fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
