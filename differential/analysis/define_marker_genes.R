library(ggpubr)

#########
## I/O ##
#########

source("/Users/ricard/gastrulation10x/settings.R")
io$diff.dir <- paste0(io$basedir,"/results/differential")
io$outdir <- paste0(io$basedir,"/results/marker_genes"); dir.create(io$outdir, showWarnings = F)

#############
## Options ##
#############

opts$groups <- opts$celltypes.1

##################
## Load results ##
##################

dt <- list()
for (i in opts$groups) {
  for (j in opts$groups) {
    if (i!=j) {
      file <- sprintf("%s/%s_vs_%s.txt.gz", io$diff.dir,i,j)
      if (file.exists(file)) {
        dt[[file]] <- fread(file) %>%
          .[,c("p.value","padj_fdr"):=NULL] %>%
          .[,c("groupA","groupB"):=list(i,j)] %>%
          setnames(sprintf("detection_rate_%s",i),"detection_rate_groupA") %>%
          setnames(sprintf("detection_rate_%s",j),"detection_rate_groupB") %>%
          .[,comparison:=paste(i,j,sep="_vs_")]
      } else {
        print(sprintf("%s does not exist",file))
      }
    }
  }
}
dt <- rbindlist(dt)

dt[is.na(sig),sig:=F]

dt[,direction:=c("down","up")[as.numeric(logFC<0)+1]]
# dt[,direction:=c("down","up")[as.numeric(detection_rate_groupA>detection_rate_groupB)+1]]

#########################
## Define marker genes ##
#########################

# Minimum fraction of significant differential pairwisecomparisons
opts$score <- 0.75

dt.filt <- dt[,.(score=round(mean(sig==T & direction=="up"),3)), by=c("groupA","gene","ens_id")] %>%
  .[score>=opts$score] %>%
  setnames("groupA","celltype") %>%
  setorder(celltype,-score)

# foo[groupA=="Mesoderm"] %>% View

##########
## Plot ##
##########

# Plot number of marker genes per cell types

to.plot <- dt.filt %>% .[,.N,by=c("celltype")]

pdf(sprintf("%s/pdf/barplot_number_marker_genes.pdf",io$outdir), width = 9, height = 5)
ggbarplot(to.plot, x="celltype", y="N", fill="celltype") +
  scale_fill_manual(values=opts$celltype.colors.1) +
  labs(x="", y="Number of marker genes") +
  theme(
    axis.text.y = element_text(size=rel(0.75)),
    axis.text.x = element_text(colour="black",size=rel(0.8), angle=90, hjust=1, vjust=0.5),
    axis.ticks.x = element_blank(),
    legend.position = "none"
)
dev.off()

# Plot exclusivity of cell types
to.plot <- dt.filt %>% .[,N:=.N,by="gene"]

pdf(sprintf("%s/pdf/boxplot_exclusivity.pdf",io$outdir), width = 9, height = 5)
ggboxplot(to.plot, x="celltype", y="N", fill="celltype", color="black") +
  scale_fill_manual(values=opts$celltype.colors.1) +
  labs(x="", y="Exclusivity of gene markers\n(the smaller the more exclusive)") +
  theme(
    axis.text.y = element_text(size=rel(0.75)),
    axis.title.y = element_text(size=rel(0.85)),
    axis.text.x = element_text(colour="black",size=rel(0.7), angle=90, hjust=1, vjust=0.5),
    legend.position = "none"
  )
dev.off()

# Plot exclusivity of marker genes
to.plot <- dt.filt %>%
  .[,.(Nx=.N),by="gene"] %>%
  .[,Nx:=factor(Nx)] %>%
  .[,.(Ny=.N),by="Nx"]

pdf(sprintf("%s/pdf/boxplot_exclusivity2.pdf",io$outdir), width = 7, height = 5)
ggbarplot(to.plot, x="Nx", y="Ny", fill="gray70") +
  labs(x="Number of different cell types per marker gene", y="") +
  theme(
    axis.text = element_text(size=rel(0.75)),
  )
dev.off()

##########
## Save ##
##########

fwrite(dt.filt, paste0(io$outdir,"/marker_genes.txt.gz"))
