source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$sce.pseudobulk <- file.path(io$basedir,"results/pseudobulk/SingleCellExperiment_pseudobulk_celltype.rds")
io$outdir <- file.path(io$basedir,"results/differential/celltypes/TFs/pseudobulk"); dir.create(io$outdir, showWarnings = F)
io$TFs <- file.path(io$basedir,"results/differential/celltypes/TFs/TFs.txt")

# Options
opts$min.cells <- 25


##############################
## Load pseudobulk RNA data ##
##############################

sce <- readRDS(io$sce.pseudobulk)[,opts$celltypes]

##################
## Rename genes ##
##################

gene_metadata <- fread(io$gene_metadata) %>% .[,c("ens_id","symbol")] %>% .[symbol!="" & ens_id%in%rownames(sce)]
sce <- sce[rownames(sce)%in%gene_metadata$ens_id,]
foo <- gene_metadata$symbol; names(foo) <- gene_metadata$ens_id
rowData(sce)$gene <- foo[rownames(sce)]
stopifnot(!is.na(rowData(sce)$gene))

################
## Subset TFs ##
################

TFs <- fread(io$TFs)[[1]] %>% str_to_title
sce <- sce[rowData(sce)$gene%in%TFs]

###############################################
## Differential expression between celltypes ##
###############################################

# i <- "Gut"; j <- "NMP"
for (i in opts$celltypes) {
  for (j in  opts$celltypes) {
    tmp <- data.table(
      ens_id = rownames(sce),
      gene = rowData(sce)$gene,
      expr_A = logcounts(sce[,i])[,1] %>% round(2),
      expr_B = logcounts(sce[,j])[,1] %>% round(2)
    ) %>% .[,diff:=round(expr_B-expr_A,2)] %>% sort.abs("diff") 

    # save      
    fwrite(tmp, sprintf("%s/%s_vs_%s.txt.gz", io$outdir,i,j), sep="\t")
  }
}
