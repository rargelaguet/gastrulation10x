
library(GGally)
library(network)
library(sna)
library(ggraph)
library(igraph)
library(tidygraph)

paga_connectivity.mtx <- fread(io$paga.connectivity) %>%
  matrix.please %>% .[opts$celltypes,opts$celltypes]

paga_coordinates.mtx <- fread(io$paga.coordinates) %>% 
  matrix.please %>% .[opts$celltypes,]

# Parse data
paga_connectivity.mtx[paga_connectivity.mtx<0.20] <- 0
paga_connectivity.mtx[paga_connectivity.mtx>=0.20] <- 1

# Create igraph object
igraph.paga <- graph_from_adjacency_matrix(paga_connectivity.mtx, mode = "undirected")

# Create tbl_graph object
igraph.paga.tbl <- as_tbl_graph(igraph.paga) %>%
  activate(nodes) %>%
  mutate(celltype=rownames(paga_connectivity.mtx)) %>%
  mutate(x=paga_coordinates.mtx[,1]) %>% mutate(y=paga_coordinates.mtx[,2])

# Create network object
net.paga = network(paga_connectivity.mtx)
net.paga %v% "x" = paga_connectivity.mtx[, 1]
net.paga %v% "y" = paga_connectivity.mtx[, 2]

##########
## TEST ##
##########

# sum(paga_connectivity.mtx==1)
# paga_connectivity.mtx["Epiblast","Rostral_neurectoderm"]
# igraph.paga.tbl %>% activate(edges) %>% as.data.table()  %>% nrow
# filter(celltype=="Epiblast")
