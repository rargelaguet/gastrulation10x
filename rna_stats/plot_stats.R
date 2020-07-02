library(ggpubr)

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation10x/settings.R")
}

io$outdir <- paste0(io$basedir,"/results/general_stats/pdf")


###############################################
## Boxplots of general statistics per sample ##
###############################################

to.plot <- sample_metadata %>% 
  melt(id.vars=c("cell","sample","stage"), measure.vars=c("nCount_RNA","nFeature_RNA"))

p <- ggboxplot(to.plot, x = "sample", y = "value", outlier.shape=NA) +
  yscale("log10", .format = TRUE) +
  labs(x="", y="") +
  facet_wrap(~variable, scales="free_y", nrow=2) +
  theme(
    legend.position = "right",
    legend.title = element_blank()
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank()
  )

pdf(paste0(io$outdir,"/general_stats_per_sample.pdf"), width=12, height=6, useDingbats = F)
print(p)
dev.off()

##################################################
## Boxplots of general statistics per cell type ##
##################################################

to.plot <- sample_metadata %>% 
  melt(id.vars=c("cell","celltype"), measure.vars=c("nCount_RNA","nFeature_RNA"))

p <- ggboxplot(to.plot, x = "celltype", y = "value", fill="celltype", outlier.shape=NA) +
  yscale("log10", .format = TRUE) +
  labs(x="", y="") +
  scale_fill_manual(values=opts$celltype.colors) +
  facet_wrap(~variable, scales="free_y") +
  guides(fill = guide_legend(override.aes = list(size=0.25), ncol=1)) +
  scale_size(guide = 'none') +
  theme(
    legend.position = "right",
    legend.text = element_text(size=rel(0.75)),
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

pdf(paste0(io$outdir,"/general_stats_per_celltype.pdf"), width=15, height=10, useDingbats = F)
print(p)
dev.off()

