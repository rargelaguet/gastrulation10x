## Option 1 ##
# dt <- fread("/Users/ricard/data/gastrulation10x/results/marker_genes/avg_expr_per_celltype_and_gene.txt.gz") %>%
#   .[,gene:=NULL] %>% unique
# m <- dt %>% 
#   .[ens_id%in%unique(marker_genes$ens_id)] %>%
#   dcast(ens_id~group,value.var="mean_expr") %>% matrix.please %>% t
# h <- hclust(dist(m))


# Load marker genes
marker_genes <- fread("/Users/ricard/data/gastrulation10x/results/marker_genes/marker_genes.txt.gz") %>%
  .[,celltype:=stringr::str_replace_all(celltype,"_", " ")] %>%
  .[celltype%in%opts$celltypes] %>%
  .[,head(.SD,n=50),by="celltype"]
table(marker_genes$celltype)