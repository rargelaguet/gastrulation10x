#########
## I/O ##
#########

source("/Users/ricard/gastrulation10x/settings.R")
io$diff.dir <- paste0(io$basedir,"/results/differential/celltypes/all_stages")
# io$outdir <- paste0(io$basedir,"/results/differential/celltypes/all_stages/pdf"); dir.create(io$outdir, showWarnings = F)

#############
## Options ##
#############

# Define cell types
# opts$groups <- opts$celltypes.1
opts$groups <- c("Caudal_epiblast","NMP")

##################
## Load results ##
##################

dt <- list()
# for (i in 1:(length(opts$groups)-1)) {
#   groupA <- opts$groups[[i]]
#   for (j in (i+1):length(opts$groups)) {
#     groupB <- opts$groups[[j]]
for (groupA in opts$groups) {
  for (groupB in opts$groups) {
    if (groupA!=groupB) {
      file <- sprintf("%s/%s_vs_%s.txt.gz", io$diff.dir,groupA,groupB)
      if (file.exists(file)) {
        dt[[file]] <- fread(file) %>%
          .[,c("groupA","groupB"):=list(groupA,groupB)] %>%
          setnames(sprintf("detection_rate_%s",groupA),"detection_rate_groupA") %>%
          setnames(sprintf("detection_rate_%s",groupB),"detection_rate_groupB") %>%
          .[,comparison:=paste(groupA,groupB,sep="_vs_")]
      } else {
        print(sprintf("%s does not exist",file))
      }
    }
  }
}
dt <- rbindlist(dt)

##########
## Plot ##
##########

total_genes <- length(unique(dt$ens_id))

# Tile plot with the fraction of differentially expressed genes between pairs of cell types
to.plot <- dt %>%
  .[,.(sum_genes=sum(sig,na.rm=T)),by=c("groupA","groupB")] %>%
  .[,fraction_genes:=sum_genes/total_genes] %>%
  .[,groupA:=factor(groupA,levels=opts$groups)] %>%
  .[,groupB:=factor(groupB,levels=opts$groups)]

p <- ggplot(to.plot, aes(x=groupA, y=groupB)) +
  geom_tile(aes(fill=fraction_genes)) +
  labs(x="", y="") +
  scale_fill_gradient(low = "white", high = "red") +
  theme_classic() +
  theme(
    axis.text.y = element_text(colour="black",size=rel(1.0)),
    axis.text.x = element_text(colour="black",size=rel(1.0), angle=90, hjust=1, vjust=0.5),
    axis.ticks = element_blank()
  )

pdf(sprintf("%s/heatmap_fraction_differential_genes.pdf",io$outdir), width = 9, height = 7)
print(p)
dev.off()

## barplots ##

to.plot <- dt %>%
  .[,.(sum_genes=sum(sig,na.rm=T)),by=c("groupA","groupB")]

for (i in unique(to.plot$groupA)) {
  p <- ggplot(to.plot[groupA==i], aes(x = reorder(groupB, -sum_genes), y = sum_genes)) + 
    facet_wrap(~groupA) +
    geom_bar(stat = "identity") +
    theme_classic() +
    labs(x="", y="Number of diff. genes") +
    theme(
      axis.text.y = element_text(colour="black",size=rel(1.0)),
      axis.text.x = element_text(colour="black",size=rel(1.0), angle=90, hjust=1, vjust=0.5),
      axis.ticks.x = element_blank()
    )
  
  pdf(sprintf("%s/%s_barplots_numbers_differential_genes.pdf",io$outdir,i), width = 7, height = 5)
  print(p)
  dev.off()
}


##########
## TEST ##
##########

to.plot %>% setorder(sum_genes) %>% head(n=25)

to.plot[groupB=="PGC" | groupA=="PGC"]

foo <- dt[groupA=="Mixed_mesoderm" & groupB=="PGC" & sig==T] %>% setorder(-log_padj_fdr)
foo <- dt[groupA=="Visceral_endoderm" & groupB=="Primitive_Streak" & sig==T] %>% setorder(-log_padj_fdr)


dt[sig==T & abs(logFC)>1.5] %>% View
