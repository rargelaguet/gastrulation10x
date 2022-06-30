#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

io$differential_results <- file.path(io$basedir,"results/differential/celltypes/TFs/pseudobulk")
io$outdir <- file.path(io$basedir,"results/differential/celltypes/TFs/pseudobulk/out"); dir.create(io$outdir, showWarnings = F)

##################
## Load results ##
##################

# i <- "Gut"; j <- "NMP"
dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- file.path(io$differential_results,sprintf("%s_vs_%s.txt.gz",i,j))
  if (file.exists(file)) {
    fread(file) %>%
      .[,c("celltypeA","celltypeB"):=list(i,j)] %>%
      setnames("diff","logFC") %>%
      return
  } }) %>% rbindlist }) %>% rbindlist %>%
  .[,celltypeA:=factor(celltypeA,levels=opts$celltypes)] %>%
  .[,celltypeB:=factor(celltypeB,levels=opts$celltypes)] %>%
  .[,direction:=as.factor(c("down","up"))[as.numeric(logFC<0)+1]]

##########################
## Create summary table ##
##########################

tmp <- dt %>%
  .[,.(logFC=round(mean(logFC),2)), by=c("celltypeB","gene","ens_id")] %>%
  # .[score>=opts$score] %>%
  setnames("celltypeB","celltype") %>%
  setorder(celltype,-logFC)

##########
## Save ##
##########

# Save marker score for all combination of genes and cell types
length(unique(tmp$gene))
length(unique(tmp$celltype))
fwrite(tmp, file.path(io$outdir,"differential_pseudobulk_logFC.txt.gz"), sep="\t", na="NA", quote=F)

