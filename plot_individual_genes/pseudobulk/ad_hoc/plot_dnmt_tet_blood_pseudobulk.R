source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$outdir <- file.path(io$basedir,"results/individual_genes/pseudobulk")
io$sce.pseudobulk <- file.path(io$basedir,"results/pseudobulk/SingleCellExperiment_pseudobulk_celltype.rds")

# Define options
opts$celltypes = c(
  "Haematoendothelial_progenitors",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3"
)
opts$min.cells <- 10

# opts$rename <- c(
#   "Haematoendothelial_progenitors" = "Haemato._progenitors"
# )

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- readRDS(io$sce.pseudobulk)[,opts$celltypes]

# Rename to gene names
gene_metadata <- fread(io$gene_metadata) %>% .[symbol!="" & ens_id%in%rownames(sce)]
sce <- sce[rownames(sce)%in%gene_metadata$ens_id,]
foo <- gene_metadata$symbol; names(foo) <- gene_metadata$ens_id
new.names <- foo[rownames(sce)]
stopifnot(sum(is.na(new.names))==0)
stopifnot(sum(duplicated(new.names))==0)
rownames(sce) <- new.names

###############
## Lineplots ##
###############

genes.to.plot <- c("Dnmt1","Dnmt3a","Dnmt3b","Tet1","Tet2","Tet3")

to.plot <- logcounts(sce)[genes.to.plot,]  %>%
  reshape2::melt() %>% as.data.table() %>%
  setnames(c("gene","celltype","expr")) %>%
  .[,gene_class:=ifelse(grepl("Dnmt",gene),"DNMTs","TETs")]

to.plot %>% .[,celltype:=stringr::str_replace_all(celltype,opts$rename)] 

ggline(to.plot, x="celltype", y="expr", group="gene") +
  # geom_bar(stat="identity", color="black") +
  scale_color_brewer(palette="Dark2") +
  facet_wrap(~gene_class, scales="free_y") +
  labs(x="",y="Gene expression") +
  theme_classic() +
  guides(x = guide_axis(angle = 90)) +
  ggrepel::geom_text_repel(aes_string(label="gene"), nudge_x=0.25, data=to.plot[celltype=="Erythroid3"]) +
  # geom_text(aes_string(label="gene"), position=position_nudge(x=0.25), data=to.plot[celltype=="late_Erythroid"]) +
  theme(
    strip.text = element_text(size=rel(0.85)),
    axis.text.x = element_text(colour="black",size=rel(0.9)),
    axis.text.y = element_text(colour="black",size=rel(0.9)),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(colour="black",size=rel(1.0)),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size=rel(0.85))
  )

# pdf(outfile, width=10, height=9)
# print(p)
# dev.off()

