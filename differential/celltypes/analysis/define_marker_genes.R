#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

io$differential_results <- file.path(io$basedir,"results/differential/celltypes/genes/all_stages")
io$outdir <- file.path(io$basedir,"results/differential/celltypes/genes/all_stages/marker_genes"); dir.create(io$outdir, showWarnings = F)

#############
## Options ##
#############

opts$celltypes <- c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)# %>% head(n=4)

# Minimum fraction of significant differential pairwise comparisons
opts$score <- 0.75

# Significance thresholds
opts$min.log2FC <- 1
opts$fdr <- 0.01

##################
## Load results ##
##################

# i <- "Gut"; j <- "NMP"
dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- file.path(io$differential_results,sprintf("%s_vs_%s.txt.gz",i,j))
  if (file.exists(file)) {
    fread(file, select = c(1,8,2,4)) %>%
      .[gene!=""] %>%
      .[is.na(logFC),c("logFC","padj_fdr"):=list(0,1)] %>%
      .[,sig:=abs(logFC)>=opts$min.log2FC & padj_fdr<=opts$fdr] %>%
      .[,c("celltypeA","celltypeB"):=list(i,j)] %>%
      return
  } }) %>% rbindlist }) %>% rbindlist %>%
  .[,celltypeA:=factor(celltypeA,levels=opts$celltypes)] %>%
  .[,celltypeB:=factor(celltypeB,levels=opts$celltypes)] %>%
  .[,direction:=as.factor(c("down","up"))[as.numeric(logFC<0)+1]]

ncelltypes <- unique(c(as.character(unique(dt$celltypeA)),as.character(unique(dt$celltypeB)))) %>% length

# Check for duplicated genes
# dt[,.N,c("gene","celltypeA","celltypeA")] %>% View

# Filter out duplicated genes
dt <- dt[,N:=length(celltypeB),by=c("gene","celltypeA")] %>% .[N==(ncelltypes-1)] %>% .[,N:=NULL]

stop("TO-DO: CHECK THAT NO CELL TYPE COMPARISONS ARE MISSING")

#########################
## Define marker genes ##
#########################

# No need for this, we ran all pair-wise comparisons
# foo <- dt[,.(score=sum(direction=="up" & sig==TRUE,na.rm=T)), by=c("celltypeA","gene","ens_id")] %>% setnames("celltypeA","celltype")
# bar <- dt[,.(score=sum(direction=="down" & sig==TRUE,na.rm=T)), by=c("celltypeB","gene","ens_id")] %>% setnames("celltypeB","celltype")
# markers_genes.dt <- merge(foo,bar,by=c("celltype","gene","ens_id"), all=TRUE) %>% 
#   .[is.na(score.x),score.x:=0] %>% .[is.na(score.y),score.y:=0] %>%
#   .[,score:=score.x+score.y] %>%
#   .[,c("score.x","score.y"):=NULL] %>%
#   .[,score:=round(score/(ncelltypes-1),2)] %>%
#   setorder(celltype,-score)
# rm(foo,bar)

markers_genes.dt <- dt %>%
  .[,.(score=round(mean(sig==T & direction=="up"),2)), by=c("celltypeA","gene","ens_id")] %>%
  # .[score>=opts$score] %>%
  setnames("celltypeA","celltype") %>%
  setorder(celltype,-score)

# markers_genes.dt <- fread(file.path(io$basedir,"results/differential/celltypes/genes/all_stages/marker_genes/marker_genes_upregulated_all.txt.gz"))

###########################################
## Add logFC values from pseudobulk data ##
###########################################

logFC_pseudobulk.dt <- file.path(io$basedir,"results/differential/celltypes/pseudobulk/out/differential_pseudobulk_logFC.txt.gz") %>% fread

markers_genes.dt <- markers_genes.dt %>% merge(logFC_pseudobulk.dt, by=c("celltype","gene","ens_id"))

##########
## Save ##
##########

# Save marker score for all combination of genes and cell types
length(unique(markers_genes.dt$gene))
length(unique(markers_genes.dt$celltype))
fwrite(markers_genes.dt, file.path(io$outdir,"marker_genes_upregulated_all.txt.gz"), sep="\t")

# Save marker score for strong markers
markers_genes_filt.dt <- markers_genes.dt %>% .[score>=opts$score & logFC>=2]
length(unique(markers_genes_filt.dt$gene))
length(unique(markers_genes_filt.dt$celltype))
fwrite(markers_genes_filt.dt, file.path(io$outdir,"marker_genes_upregulated_filtered.txt.gz"), sep="\t")


########################################################
## Plot score vs number of marker genes per cell type ##
########################################################

seq.ranges <- seq(0.10,1,by=0.10)
names(seq.ranges) <- as.character(1:length(seq.ranges))

to.plot <- seq.ranges %>% map(function(i) {
  markers_genes.dt %>% .[score>=i] %>% .[,.N,by=c("celltype")] %>% .[,min_score:=i] %>% return
}) %>% rbindlist

p <- ggline(to.plot, x="min_score", y="N", color="celltype", plot_type = "l") +
  scale_color_manual(values=opts$celltype.colors) +
  labs(x="Minimum marker score", y="Number of marker genes") +
  theme(
    axis.text = element_text(size=rel(0.65)),
    legend.position = "none"
  )

pdf(file.path(io$outdir,"lineplot_number_marker_genes.pdf"), width = 7, height = 4)
print(p)
dev.off()

################################################
## Plot number of marker genes per cell types ##
################################################

to.plot <- markers_genes_filt.dt %>% .[,.N,by=c("celltype")]

p <- ggbarplot(to.plot, x="celltype", y="N", fill="celltype") +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="Number of marker genes") +
  theme(
    axis.text.y = element_text(size=rel(0.65)),
    axis.text.x = element_text(colour="black",size=rel(0.7), angle=90, hjust=1, vjust=0.5),
    axis.title = element_text(colour="black",size=rel(0.75)),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

pdf(file.path(io$outdir,"barplot_number_marker_genes.pdf"), width = 8, height = 4)
print(p)
dev.off()

##################################
## Plot gene marker exclusivity ##
##################################

to.plot <- markers_genes_filt.dt %>%
  .[,.(Nx=.N),by="gene"] %>%
  .[,Nx:=factor(Nx)] %>%
  .[,.(Ny=.N),by="Nx"]

p <- ggbarplot(to.plot, x="Nx", y="Ny", fill="gray70") +
  labs(x="Number of different cell types per marker TF", y="") +
  theme(
    axis.text = element_text(size=rel(0.75)),
  )
pdf(file.path(io$outdir,"boxplot_exclusivity_per_TF.pdf"), width = 7, height = 5)
print(p)
dev.off()

################################################
## Plot gene marker exclusivity per cell type ##
################################################

to.plot <- markers_genes_filt.dt %>% .[,N:=.N,by="gene"]

p <- ggboxplot(to.plot, x="celltype", y="N", fill="celltype", color="black") +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="Exclusivity of TF markers\n(the smaller the more exclusive)") +
  theme(
    axis.text.y = element_text(size=rel(0.75)),
    axis.title.y = element_text(size=rel(0.85)),
    axis.text.x = element_text(colour="black",size=rel(0.7), angle=90, hjust=1, vjust=0.5),
    legend.position = "none"
  )

pdf(file.path(io$outdir,"boxplot_exclusivity_per_celltype.pdf"), width = 9, height = 5)
print(p)
dev.off()


