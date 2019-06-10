library(purrr)
library(data.table)
source("/Users/ricard/data/gastrulation10x/scripts/core_functions.R")
sce <- load_data(normalise = TRUE, hvg = TRUE, corrected = FALSE, remove_doublets = FALSE)


cData <- meta %>% tibble::remove_rownames() %>% tibble::column_to_rownames("cell") %>% DataFrame
colData(sce_all) <- cData

# sce_filt <- sce_all[unlist(hvg.list),sce_all$stage=="E6.5"]

###########

# blanca_meta <- fread("/Users/ricard/data/gastrulation10x/blanca_meta.tab")

# E6.5
sce_filt <- list()
meta_filt <- list()
sce_filt[["E6.5"]] <- sce_all[hvg.list$E6.5,sce_all$stage=="E6.5"]
meta_filt[["E6.5"]] <- blanca_meta[stage=="E6.5"]
sce_filt[["E6.5"]] <- sce_filt[["E6.5"]][,meta_filt[["E6.5"]]$cell]

colData(sce_filt[["E6.5"]]) <- meta_filt[["E6.5"]] %>% tibble::remove_rownames() %>% tibble::column_to_rownames("cell") %>% DataFrame
colnames(sce_filt[["E6.5"]]) <- meta_filt[["E6.5"]]$cell

plotTSNE(sce_filt[["E6.5"]], colour_by="cluster.ann", ncomponents=2)

head(cData)
