library(ggpubr)

#########
## I/O ##
#########

source("/Users/ricard/gastrulation10x/settings.R")
io$diff.dir <- paste0(io$basedir,"/results/differential")
io$outdir <- paste0(io$basedir,"/results/cellassign"); dir.create(io$outdir, showWarnings = F)

#############
## Options ##
#############

opts <- list()

# Define cell types
# opts$groups <- c(
#   "Epiblast/PS",
#   "ExE_ectoderm",
#   "ExE_endoderm",
#   "Nascent_mesoderm",
#   "Neuroectoderm",
#   "Blood_progenitors",
#   "Mesoderm",
#   "ExE_mesoderm",
#   "Caudal_epiblast",
#   "PGC",
#   "Mesenchyme",
#   "Surface_ectoderm",
#   "Endoderm",
#   "Notochord",
#   "Erythroid",
#   "Parietal_endoderm",
#   "Endothelium",
#   "Spinal_cord",
#   "Cardiomyocytes",
#   "NMP"
# )

##################
## Load results ##
##################

dt <- list()
for (i in opts$groups) {
  i <- opts$groups[[i]]
  for (j in opts$groups) {
    if (i!=j) {
      j <- opts$groups[[j]]
      file <- sprintf("%s/%s_vs_%s.txt.gz", io$diff.dir,i,j)
      if (file.exists(file)) {
        dt[[file]] <- fread(file) %>%
          .[,c("p.value","padj_fdr"):=NULL] %>%
          .[,c("groupA","groupB"):=list(i,j)] %>%
          setnames(sprintf("detection_rate_%s",i),"detection_rate_groupA") %>%
          setnames(sprintf("detection_rate_%s",j),"detection_rate_groupB") %>%
          .[,comparison:=paste(i,j,sep="_vs_")]
      } else {
        sprintf("%s does not exist",file)
      }
    }
  }
}
dt <- rbindlist(dt)

dt[is.na(sig),sig:=F]

# dt[,direction:=c("-","+")[as.numeric(logFC<0)+1]]
dt[,direction:=c("down","up")[as.numeric(detection_rate_groupA>detection_rate_groupB)+1]]

#########################
## Define marker genes ##
#########################

# Minimum fraction of significant differential pairwisecomparisons
opts$score <- 0.75

dt.filt <- dt[,.(score=mean(sig==T & direction=="up")),by=c("groupA","gene")] %>%
  .[score>=opts$foo] %>%
  setnames("groupA","celltype") %>%
  setorder(celltype,-score)

# foo[groupA=="Mesoderm"] %>% View

##########
## Plot ##
##########

# Plot number of marker genes per cell types

to.plot <- dt.filt %>% .[,.N,by=c("celltype")]

ggbarplot(to.plot, x="celltype", y="N") +
  labs(x="", y="Number of marker genes") +
  theme(
    axis.text.y = element_text(size=rel(0.75)),
    axis.text.x = element_text(colour="black",size=rel(0.8), angle=90, hjust=1, vjust=0.5),
    axis.ticks.x = element_blank()
  )

# Plot exclusivty of cell types


# Plot exclusivity of marker genes

# to.plot <- dt.filt %>%
#   .[,.N,by="gene"] %>%
#   .[,as]
#   setorder(-N)
# 
# gghistogram(to.plot, x="N") +
#   labs(x="Number of cell types", y="") +
#   theme(
#     axis.text = element_text(size=rel(0.75)),
#   )

##########
## Save ##
##########

fwrite(dt.filt, paste0(io$outdir,"/marker_genes.txt.gz"))