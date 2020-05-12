matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}

dt <- fread("/Users/ricard/data/gastrulation10x/results/marker_genes/avg_expr_per_celltype_and_gene.txt.gz") %>%
  .[,gene:=NULL] %>% unique

marker_genes <- fread("/Users/ricard/data/gastrulation10x/results/marker_genes/marker_genes.txt.gz") %>%
  .[,head(.SD,n=50),by="celltype"]

m <- dt %>% 
  .[ens_id%in%unique(marker_genes$ens_id)] %>%
  dcast(ens_id~group,value.var="mean_expr") %>% matrix.please %>% t

dim(m)
h <- hclust(dist(m))

plot(h)

k <- 2
cutree(h,k=k)


##########
## Plot ##
##########


hcd <- as.dendrogram(h)
plot(hcd, type = "rectangle", ylab = "Height")
abline(h=h$height[length(h$height)-(k-1)],col="red")


library("ggplot2")
library("ggdendro")

ggdendrogram(h)
