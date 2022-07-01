
#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

io$outdir <- file.path(io$basedir,"results/gene_statistics")

###############
## Load data ##
###############

# Load cell metadata
sample_metadata <- fread(args$metadata) %>%
  .[stripped==F & doublet==F & !is.na(celltype)]

# Load SingleCellExperiment, single cells
sce_cells <- load_SingleCellExperiment(io$sce, cells=sample_metadata$cell, normalise = T)
colData(sce_cells) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

# Load SingleCellExperiment, pseudobulk by cell type
sce_pseudobulk <- readRDS(io$sce_pseudobulk_celltype)

##################
## Calculations ##
##################

stopifnot(rownames(sce_cells)==rownames(sce_pseudobulk))

dt <- data.table(
  ens_id = rownames(sce_cells),
  mean_single_cells = Matrix::rowMeans(logcounts(sce_cells)) %>% round(2),
  var_single_cells = sparseMatrixStats::rowVars(logcounts(sce_cells)) %>% round(2),
  mean_pseudobulk = apply(logcounts(sce_pseudobulk),1,mean) %>% round(2),
  var_pseudobulk = apply(logcounts(sce_pseudobulk),1,var) %>% round(2)
) %>% merge(fread(io$gene_metadata)[,c("ens_id","symbol")] %>% setnames("symbol","gene"), by="ens_id",all.y=T)

##########
## Plot ##
##########

ggscatter(dt, x="mean_single_cells", y="mean_pseudobulk", size=0.5)
ggscatter(dt, x="var_single_cells", y="var_pseudobulk", size=0.5)


ggscatter(dt, x="mean_single_cells", y="var_single_cells", size=0.5)
ggscatter(dt, x="mean_pseudobulk", y="var_pseudobulk", size=0.5)


#############
## Explore ##
#############

i <- "ENSMUSG00000022484"

dt[var_pseudobulk<2 & var_single_cells>2] %>% View
to.plot <- logcounts(sce_pseudobulk[i,])[1,] %>% as.data.table(keep.rownames = T) %>%
  setnames(c("celltype","x"))

sort(counts(sce_pseudobulk[i,])[1,])["Epiblast"]
sort(counts(sce_pseudobulk[i,])[1,])["Notochord"]
sort(logcounts(sce_pseudobulk[i,])[1,])["Epiblast"]
sort(logcounts(sce_pseudobulk[i,])[1,])["Notochord"]

ggbarplot(to.plot, x="celltype", y="x", fill="celltype") +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="Expression") +
  guides(x = guide_axis(angle = 90)) +
  theme(
    legend.position = "none",
    # axis.text.x = element_blank(),
    axis.text.x = element_text(size=rel(0.75)),
    axis.text.y = element_text(size=rel(0.75)),
    axis.ticks.x = element_blank()
  )

# to.plot <- logcounts(sce_cells["ENSMUSG00000012396",])[1,] %>% as.data.table(keep.rownames = T) %>%
#   setnames(c("celltype","x"))
# 
# ggboxplot(to.plot, x="celltype", y="x", fill="celltype") +
#   scale_fill_manual(values=opts$celltype.colors) +
#   labs(x="", y="Expression") +
#   guides(x = guide_axis(angle = 90)) +
#   theme(
#     legend.position = "none",
#     # axis.text.x = element_blank(),
#     axis.text.x = element_text(size=rel(0.75)),
#     axis.text.y = element_text(size=rel(0.75)),
#     axis.ticks.x = element_blank()
#   )

##########
## Save ##
##########

fwrite(dt, file.path(io$outdir,"gene_statistics.txt.gz"))
