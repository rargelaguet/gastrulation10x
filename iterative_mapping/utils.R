
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
    # print(table(sce.sub$celltype))
    
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
  groups.parsed <- c(paste(groupA,collapse="%"), paste(groupB,collapse="%"))
  sce.atlas$group <- groups.parsed[as.numeric(sce.atlas$celltype%in%groupA)+1]
  sce.atlas$group <- factor(sce.atlas$group,levels=groups.parsed)
  
  # Run differential expression between groupA and groupB cells
  diff <- differential_expression(sce.atlas, groupA, groupB)
  
  # Select gene markers for the two group
  markers.groupA <- diff[padj_fdr<opts$threshold.fdr & logFC>opts$min.logFC,ens_id] %>% head(n=opts$number.diff.genes)
  markers.groupB <- diff[padj_fdr<opts$threshold.fdr & logFC<(-opts$min.logFC),ens_id] %>% head(n=opts$number.diff.genes)
  marker_list <- list(markers.groupA,markers.groupB)
  names(marker_list) <- groups.parsed
  
  # Run cellassign
  cellassign.fit <- cellassign.fn(sce.query, marker_list)
  
  # Fetch cell type predictions
  pred.celltypes <- celltypes(cellassign.fit, assign_prob = 0.50)
  
  # Return data.table
  dt <- data.table(cell=colnames(sce.query), celltype.pred=pred.celltypes)
  
  return(dt)
    
}



differential_expression <- function(sce.atlas, groupA, groupB) {
  
  
  # Filter genes by detection rate per group
  cdr_A <- rowMeans(logcounts(sce.atlas[,sce.atlas$celltype%in%groupA])>0) >= opts$min_detection_rate_per_group
  cdr_B <- rowMeans(logcounts(sce.atlas[,sce.atlas$celltype%in%groupB])>0) >= opts$min_detection_rate_per_group
  sce.atlas <- sce.atlas[cdr_A|cdr_B,]
  
  # Convert SCE to DGEList
  sce_edger <- scran::convertTo(sce.atlas, type="edgeR")
  
  # Define design matrix
  cdr <- colMeans(logcounts(sce.atlas)>0)
  design <- model.matrix(~cdr+sce.atlas$group)
  
  # Estimate dispersions
  sce_edger <- estimateDisp(sce_edger,design)
  
  # Fit GLM
  fit <- glmQLFit(sce_edger,design)
  
  # Likelihood ratio test
  lrt <- glmQLFTest(fit)
  
  # Construct output data.frame
  out <- topTags(lrt, n=nrow(lrt))$table %>% 
    as.data.table(keep.rownames=T) %>%
    setnames(c("ens_id","logFC","logCPM","LR","p.value","padj_fdr")) %>%
    .[,log_padj_fdr:= -log10(padj_fdr)] %>%
    .[,c("logCPM","LR"):=NULL] %>%
    setorder(padj_fdr)
  
  return(out)
}

cellassign.fn <- function(sce, marker_list) {
  
  # Create binary membership matrix  
  bmat <- marker_list_to_mat(marker_list)
  bmat <- bmat[,colnames(bmat)!="other"]
  
  # Subset genes
  sce <- sce[rownames(bmat),]
  
  # Extract size factors
  s <- sizeFactors(sce)
  
  # Run
  fit <- cellassign(sce, 
    marker_gene_info = bmat, 
    min_delta = 1,
    s = s, 
    # learning_rate = 1e-2, 
    shrinkage = TRUE,
    verbose = FALSE
  )
  return(fit)
}
