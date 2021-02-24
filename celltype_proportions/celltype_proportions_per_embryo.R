#####################
## Define settings ##
#####################

source("/Users/ricard/gastrulation10x/settings.R")
io$outdir <- paste0(io$basedir,"/results/celltype_proportions")

################
## Parse data ##
################

sample_metadata[,.N,by=c("stage","sample")] %>% setorder(-stage) %>% head

# Remove some lineages
# sample_metadata <- sample_metadata[!celltype%in%c("ExE_ectoderm","ExE_endoderm","Parietal_endoderm","Visceral_endoderm")]

##################
## Calculations ##
##################

dt <- sample_metadata %>% 
  .[,N:=.N,by="sample"] %>%
  .[,.(N=.N, celltype_proportion=.N/unique(N)),by=c("sample","stage","celltype")] %>%
  setorder(sample)

##########
## Save ##
##########

# fwrite(dt, paste0(io$outdir,"/celltype_proportions.txt.gz"), sep="\t")

##########
## Plot ##
##########

celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% dt$celltype]
stopifnot(sort(unique(as.character(dt$celltype))) == sort(names(celltype.colors)))

to.plot <- dt %>% 
  .[,stage:=gsub("\\.","_",stage)] %>%
  .[, celltype:=factor(celltype,levels=sort(names(celltype.colors), decreasing = F))]

for (i in unique(dt$stage)) {
  
  p <- ggplot(to.plot[stage==i], aes(x=celltype, y=N)) +
    geom_bar(aes(fill=celltype), stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors, drop=FALSE) +
    scale_x_discrete(drop=FALSE) +
    facet_wrap(~sample, nrow=1, scales="fixed") +
    coord_flip() +
    labs(y="Number of cells") +
    theme_classic() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(1.0)),
      axis.title.x = element_text(color="black", size=rel(1.0)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1.0), color="black"),
      axis.text.x = element_text(size=rel(1.0), color="black")
    )
  
  pdf(sprintf("%s/pdf/%s.pdf",io$outdir,i))
  # pdf(sprintf("%s/pdf/%s.pdf",io$outdir,i), width=10, height=5.5, useDingbats = F)
  print(p)
  dev.off()
}
  