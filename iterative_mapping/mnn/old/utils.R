
recursive.fn <- function(sce.query, sce.atlas, dist) {
  
  cell.types.to.loop <- sce.query$pred.celltype %>% unique %>% .[grep("%",.)]
  
  dt_list <- list()
  for (i in cell.types.to.loop) {
    ids <- which(sce.query$pred.celltype==i)
    foo <- stringr::str_split(i,"%") %>% unlist
    
    # Subset
    dist.sub <- usedist::dist_subset(dist, foo)
    h.sub <- hclust(dist.sub, method="ward.D")
    sce.query.sub <- sce.query[,ids]
    sce.atlas.sub <- sce.atlas[,sce.atlas$celltype %in% foo]
    
    # Run cell type assignment prediction
    dt_list[[i]] <- mapping.fn(sce.query.sub, sce.atlas.sub, h.sub)
  }
  dt <- rbindlist(dt_list)
  return(dt)
    
}

mapping.fn <- function(sce.query, sce.atlas, h) {
  
  # Cut the tree into two groups
  cut <- cutree(h, k=2)
  
  # Define groups
  groupA <- names(which(cut==1))
  groupB <- names(which(cut==2))
  groups.parsed <- c(paste(groupB,collapse="%"),paste(groupA,collapse="%"))
  sce.atlas$group <- groups.parsed[as.numeric(sce.atlas$celltype%in%groupA)+1]
  # sce.atlas$group <- factor(sce.atlas$group,levels=groups.parsed)
  
  # Run MNN
  dt <- mnn.fn(sce.query, sce.atlas)
  
  return(dt)
    
}


mnn.fn <- function(sce.query, sce.atlas, npcs = 25, k = 15) {

  # Sanity checks

  # Concatenate
  sce_all <- SingleCellExperiment(
    list(counts=Matrix::Matrix(cbind(counts(sce.atlas),counts(sce.query)), sparse=TRUE))
  )
  block <- c(rep("atlas",ncol(sce.atlas)),rep("query",ncol(sce.query)))

  # Normalize
  sce_all <- multiBatchNorm(sce_all, batch=as.factor(block))
  # sce_all <- logNormCounts(sce_all)

  # Computing highly variable genes
  hvgs <- getHVGs(sce_all, block = block, min.mean = 1e-3, p.value = 0.01)

  #  PCA
  # pca_all <- irlba::prcomp_irlba(t(logcounts(sce_all)), n = npcs)$x
  pca_all <- multiBatchPCA(sce_all,
    batch = block,
    subset.row = hvgs,
    d = npcs,
    preserve.single = TRUE,
    assay.type = "logcounts"
  )[[1]]
  rownames(pca_all) <- colnames(sce_all)
  atlas_pca <- pca_all[1:ncol(sce.atlas),]
  query_pca <- pca_all[-(1:ncol(sce.atlas)),]


  # MNN mapping
  # correct <- reducedMNN(rbind(atlas_pca, query_pca),
  correct <- reducedMNN(pca_all, batch = block)[["corrected"]]
  correct_atlas <- correct[1:nrow(atlas_pca),]
  correct_query   <- correct[-(1:nrow(atlas_pca)),]

  # get metadata
  meta_query <- as.data.frame(colData(sce.query)) %>% tibble::rownames_to_column("cell") %>% .[,c("cell"),drop=F]
  meta_atlas <- as.data.frame(colData(sce.atlas)) %>% tibble::rownames_to_column("cell") %>% .[,c("cell","group")]
  mapping <- get_meta(
    correct_atlas = correct_atlas,
    atlas_meta = meta_atlas,
    correct_query = correct_query,
    query_meta = meta_query,
    k_map = k
  )


  # Computing mapping scores
  # out <- list()
  # for (i in seq(from = 1, to = k)) {
  #   out$closest.cells[[i]]     <- sapply(mapping, function(x) x$cells.mapped[i])
  #   out$celltypes.mapped[[i]]  <- sapply(mapping, function(x) x$celltypes.mapped[i])
  # }  
  # multinomial.prob <- getMappingScore(out)

  # ## 6. Prepare output
  # out$correct_atlas <- correct_atlas
  # out$correct_query <- correct_query
  ct <- sapply(mapping, function(x) x$celltype.mapped); is.na(ct) <- lengths(ct) == 0
  cm <- sapply(mapping, function(x) x$cells.mapped[1]); is.na(cm) <- lengths(cm) == 0

  # Return data.table
  mapping.dt <- data.table(
    cell = names(mapping), 
    celltype.pred = unlist(ct)
    # closest.cell    = unlist(cm)
    # score = multinomial.prob[[1]]
  )  
  return(mapping.dt)
}


get_meta <- function(correct_atlas, atlas_meta, correct_query, query_meta, k_map = 10){
  knns <- BiocNeighbors::queryKNN(correct_atlas, correct_query, k = k_map, get.index = TRUE, get.distance = FALSE)
  
  #get closest k matching cells
  k.mapped  <- t(apply(knns$index, 1, function(x) atlas_meta$cell[x]))
  
  # get celltypes
  celltypes <- t(apply(k.mapped, 1, function(x) atlas_meta$group[match(x, atlas_meta$cell)]))
  celltype.mapped <- apply(celltypes, 1, function(x) getmode(x, 1:length(x)))
  
  out <- lapply(1:length(celltype.mapped), function(x){
    list(cells.mapped     = k.mapped[x,],
         celltype.mapped  = celltype.mapped[x],
         celltypes.mapped = celltypes[x,]
        )
  })
  names(out) <- query_meta$cell
  return(out)
}


getmode <- function(v, dist) {
  tab <- table(v)
  #if tie, break to shortest distance
  if(sum(tab == max(tab)) > 1){
    tied <- names(tab)[tab == max(tab)]
    sub  <- dist[v %in% tied]
    names(sub) <- v[v %in% tied]
    return(names(sub)[which.min(sub)])
  } else {
    return(names(tab)[which.max(tab)])
  }
}

getHVGs <- function(sce, block, min.mean = 1e-3, p.value=0.01){
  decomp <- modelGeneVar(sce, block=block)
  decomp <- decomp[decomp$mean > min.mean,]
  decomp$FDR <- p.adjust(decomp$p.value, method = "fdr")
  return(rownames(decomp)[decomp$p.value < p.value])
}
