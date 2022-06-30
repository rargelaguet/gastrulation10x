#####################
## Define settings ##
#####################

source("/Users/argelagr/gastrulation10x/settings.R")
io$outdir <- paste0(io$basedir,"/results/celltype_proportions")

opts$remove_ExE_cells <- TRUE

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>%
  .[stripped==F & doublet==F] %>%
  .[,celltype:=factor(celltype, levels=names(opts$celltype.colors))]

if (opts$remove_ExE_cells) {
  sample_metadata <- sample_metadata %>%
    .[!celltype%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

# sample_metadata[,.N,by=c("stage","sample")] %>% setorder(-stage) %>% head

##################
## Calculations ##
##################

dt <- sample_metadata %>% 
  .[,N:=.N,by="sample"] %>%
  .[,.(N=.N, celltype_proportion=round(.N/unique(N),2)),by=c("sample","stage","celltype")] %>%
  setorder(sample)

# Save
fwrite(dt, file.path(io$outdir,"celltype_proportions_noExE.txt.gz"), sep="\t")

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
  