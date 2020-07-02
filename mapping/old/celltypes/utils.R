joint.normalisation <- function(sce.query, sce.atlas, cosineNorm = FALSE) {
  
  block <- c(sce.atlas$sample,sce.query$sample) %>% as.factor
  
  if (isTRUE(cosineNorm)) {
    assay <- "cosineNorm"
    
    # Log normalisation per data set
    if (length(unique(sce.atlas$sample))>1) {
      sce.atlas <- multiBatchNorm(sce.atlas, batch=as.factor(sce.atlas$sample))
    } else {
      sce.atlas <- logNormCounts(sce.atlas)
    }
    
    if (length(unique(sce.query$sample))>1) {
      sce.query <- multiBatchNorm(sce.query, batch=as.factor(sce.query$sample))
    } else {
      sce.query <- logNormCounts(sce.query)
    }
    
    # Cosine normalisation
    assay(sce.atlas, assay) <- cosineNorm(assay(sce.atlas, "logcounts"))
    assay(sce.query, assay) <- cosineNorm(assay(sce.query, "logcounts"))
    
    # Concatenate
    sce.all <- SingleCellExperiment(
      list(cosineNorm=cbind(assay(sce.atlas,assay), assay(sce.query,assay)))
    )
    
  } else {
    assay <- "logcounts"
    
    # Concatenate
    sce.all <- SingleCellExperiment(
      list(counts=Matrix::Matrix(cbind(counts(sce.atlas),counts(sce.query)),sparse=TRUE))
    )
    
    # Log Normalise
    sce.all <- multiBatchNorm(sce.all, batch=block)
  }
  
  # Create block vector
  sce.all$block1 <- c(sce.atlas$sample,sce.query$sample) %>% as.factor
  sce.all$block <- c(rep("atlas",ncol(sce.atlas)),rep("query",ncol(sce.query))) %>% as.factor
  
  return(sce.all)
}

mnn.fn <- function(sce.all, sce.query, sce.atlas, genes = NULL, npcs = 50, k = 25, cosineNorm = FALSE) {
  
  # Note that normalisation must be done beforehand using the joint.normalisation function
  block <- sce.all$block
  
  if ("cosineNorm" %in% names(assays(sce.all))) {
    assay <- "cosineNorm"
  } else {
    assay <- "logcounts"
  }
  
  # Highly variable genes
  # (Q) Not sure if selecting HVGs using only the atlas is a good idea
  # if (is.null(genes)) {
  #   # genes <- getHVGs(sce.atlas, block=as.factor(sce.atlas$sample), assay.type="logcounts")
  #   genes <- getHVGs(sce.all, block=block, assay.type="logcounts")
  #   # genes <- rownames(sce.atlas)
  # }
  
  # PCA
  pca_all <- multiBatchPCA(sce.all,
    batch = block,
    subset.row = genes,
    d = npcs,
    preserve.single = TRUE,
    assay.type = assay
  )[[1]]
  rownames(pca_all) <- colnames(sce.all)
  atlas_pca <- pca_all[1:ncol(sce.atlas),]
  query_pca <- pca_all[-(1:ncol(sce.atlas)),]

  # (Optional) Batch effect correction for the atlas
  atlas_meta <- colData(sce.atlas)[,c("stage","sample")] %>% as.data.frame
  stages <- sort(unique(sce.atlas$stage))
  order_df        <- atlas_meta[!duplicated(atlas_meta$sample), c("stage", "sample")]
  order_df$ncells <- sapply(order_df$sample, function(x) sum(atlas_meta$sample == x))
  order_df$stage  <- factor(order_df$stage, levels = stages)
  order_df       <- order_df[order(order_df$stage, order_df$ncells, decreasing = TRUE),]
  order_df$stage <- as.character(order_df$stage)
  set.seed(42)
  atlas_corrected <- doBatchCorrect(
    sce.atlas       = sce.atlas,
    timepoints      = atlas_meta$stage,
    samples         = atlas_meta$sample,
    timepoint_order = order_df$stage,
    sample_order    = order_df$sample,
    pca             = atlas_pca
  )
  
  
  # MNN mapping
  correct <- reducedMNN(pca_all, batch = block)[["corrected"]]
  correct_atlas <- correct[1:nrow(atlas_pca),,drop=F]
  correct_query   <- correct[-(1:nrow(atlas_pca)),,drop=F]

  # get metadata
  meta_query <- as.data.frame(colData(sce.query)) %>% tibble::rownames_to_column("cell") %>% .[,c("cell"),drop=F]
  meta_atlas <- as.data.frame(colData(sce.atlas)) %>% tibble::rownames_to_column("cell") %>% .[,c("cell","celltype","stage")]
  mapping <- get_meta(
    correct_atlas = correct_atlas,
    atlas_meta = meta_atlas,
    correct_query = correct_query,
    query_meta = meta_query,
    k_map = k
  )
  
  # Compute mapping scores
  out <- list()
  for (i in seq(from = 1, to = k)) {
    out$closest.cells[[i]]     <- sapply(mapping, function(x) x$cells.mapped[i])
    out$celltypes.mapped[[i]]  <- sapply(mapping, function(x) x$celltypes.mapped[i])
    out$stages.mapped[[i]] <- sapply(mapping, function(x) x$stages.mapped[i])
  }  
  foo <- getMappingScore(out)
  celltype.multinomial.prob <- foo$celltype.score
  stage.multinomial.prob <- foo$stage.score

  #  Prepare output
  ct <- lapply(mapping, function(x) x$celltype.mapped); is.na(ct) <- lengths(ct) == 0
  st <- lapply(mapping, function(x) x$stage.mapped); is.na(st) <- lengths(st) == 0
  cm <- lapply(mapping, function(x) x$cells.mapped[1]); is.na(cm) <- lengths(cm) == 0
  mapping.dt <- data.table(
    cell = names(mapping), 
    celltype_mapped = unlist(ct),
    stage_mapped    = unlist(st),
    celltype_score = celltype.multinomial.prob,
    stage_score = stage.multinomial.prob
  )  
  
  return(mapping.dt)
}


get_meta <- function(correct_atlas, atlas_meta, correct_query, query_meta, k_map = 10){
  knns <- BiocNeighbors::queryKNN(correct_atlas, correct_query, k = k_map, get.index = TRUE, get.distance = FALSE)
  
  #get closest k matching cells
  k.mapped  <- t(apply(knns$index, 1, function(x) atlas_meta$cell[x]))
  
  # get celltypes
  celltypes <- t(apply(k.mapped, 1, function(x) atlas_meta$celltype[match(x, atlas_meta$cell)]))
  celltype.mapped <- apply(celltypes, 1, function(x) getmode(x, 1:length(x)))
  
  # get stages
  stages    <- t(apply(k.mapped, 1, function(x) atlas_meta$stage[match(x, atlas_meta$cell)]))
  stage.mapped    <- apply(stages, 1, function(x) getmode(x, 1:length(x)))

  out <- lapply(1:length(celltype.mapped), function(x){
    list(cells.mapped     = k.mapped[x,],
         celltype.mapped  = celltype.mapped[x],
         celltypes.mapped = celltypes[x,],
         stage.mapped     = stage.mapped[x],
         stages.mapped    = stages[x,]
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

getHVGs <- function(sce, block, min.mean = 1e-3, p.value = 0.01, ...){
  decomp <- modelGeneVar(sce, block=block, ...)
  decomp <- decomp[decomp$mean > min.mean,]
  decomp$FDR <- p.adjust(decomp$p.value, method = "fdr")
  return(rownames(decomp)[decomp$p.value < p.value])
}



getMappingScore <- function(mapping){
  celltypes_accrossK <- matrix(unlist(mapping$celltypes.mapped),
                               nrow=length(mapping$celltypes.mapped[[1]]),
                               ncol=length(mapping$celltypes.mapped))
  stages_accrossK <- matrix(unlist(mapping$stages.mapped),
                                nrow=length(mapping$stages.mapped[[1]]),
                                ncol=length(mapping$stages.mapped))

  # get celltype score
  celltype.score <- c()
  for (i in 1:nrow(celltypes_accrossK)) {
    p <- max(table(celltypes_accrossK[i,]))
    index <- which(table(celltypes_accrossK[i,]) == p)
    p <- p/length(mapping$celltypes.mapped)
    celltype.score <- c(celltype.score,p)
  }

  # get stage score
  stage.score <- NULL
  for (i in 1:nrow(stages_accrossK)) {
      p <- max(table(stages_accrossK[i,]))
      index <- which(table(stages_accrossK[i,]) == p)
      p <- p/length(mapping$stages.mapped)
      stage.score <- c(stage.score,p)
  }

  return(list("celltype.score"=celltype.score, "stage.score"=stage.score))  
}


doBatchCorrect <- function(sce.atlas, timepoints, samples, timepoint_order, sample_order, pca, BPPARAM = SerialParam()){
  require(BiocParallel)
  
  if(length(unique(samples)) == 1){
    return(pca)
  }
  
  #create nested list
  pc_list    <- lapply(unique(timepoints), function(tp){
    sub_pc   <- pca[timepoints == tp, , drop = FALSE]
    sub_samp <- samples[timepoints == tp]
    list     <- lapply(unique(sub_samp), function(samp){
      sub_pc[sub_samp == samp, , drop = FALSE]
    })
    names(list) <- unique(sub_samp)
    return(list)
  })
  
  names(pc_list) <- unique(timepoints)
  
  #arrange to match timepoint order
  pc_list <- pc_list[order(match(names(pc_list), timepoint_order))]
  pc_list <- lapply(pc_list, function(x){
    x[order(match(names(x), sample_order))]
  })
  
  #perform corrections within list elements (i.e. within stages)
  correct_list <- lapply(pc_list, function(x){
    if(length(x) > 1){
      return(do.call(reducedMNN, c(x, BPPARAM = BPPARAM))$corrected)
    } else {
      return(x[[1]])
    }
  })
  
  #perform correction over list
  if(length(correct_list)>1){
    correct <- do.call(reducedMNN, c(correct_list, BPPARAM = BPPARAM))$corrected
  } else {
    correct <- correct_list[[1]]
  }
  
  correct <- correct[match(colnames(sce.atlas), rownames(correct)),]
  
  return(correct)
  
}