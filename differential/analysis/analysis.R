#########
## I/O ##
#########

source("/Users/ricard/gastrulation10x/settings.R")
io$diff.dir <- paste0(io$basedir,"/results/differential")
io$outdir <- paste0(io$basedir,"/results/differential/pdf"); dir.create(io$outdir, showWarnings = F)

#############
## Options ##
#############

opts <- list()

# Define cell types
opts$groups <- c(
  "Epiblast",
  "Primitive_Streak",
  "ExE_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "Nascent_mesoderm",
  "Neurectoderm",
  "Blood_progenitors",
  "Mixed_mesoderm",
  "ExE_mesoderm",
  "Pharyngeal_mesoderm",
  "Caudal_epiblast",
  "PGC",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Surface_ectoderm",
  "Gut",
  "Paraxial_mesoderm",
  "Notochord",
  "Somitic_mesoderm",
  "Caudal_Mesoderm",
  "Erythroid",
  "Def._endoderm",
  "Parietal_endoderm",
  "Allantois",
  "Anterior_Primitive_Streak",
  "Endothelium",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Cardiomyocytes",
  "NMP",
  "Neural_crest"
)

##################
## Load results ##
##################

dt <- list()
for (i in 1:(length(opts$groups)-1)) {
  groupA <- opts$groups[[i]]
  for (j in (i+1):length(opts$groups)) {
    groupB <- opts$groups[[j]]
    file <- sprintf("%s/%s_vs-%s.txt.gz", io$diff.dir,groupA,groupB)
    if (file.exists(file)) {
      dt[[file]] <- fread(file) %>%
        .[,c("groupA","groupB"):=list(groupA,groupB)] %>%
        setnames(sprintf("detection_rate_%s",groupA),"detection_rate_groupA") %>%
        setnames(sprintf("detection_rate_%s",groupB),"detection_rate_groupB") %>%
        .[,comparison:=paste(groupA,groupB,sep="_vs_")]
    } else {
      sprintf("%s does not exist",file)
    }
  }
}
dt <- rbindlist(dt)

##########
## Plot ##
##########

total_genes <- length(unique(dt$ens_id))

# Heatmap with the fraction of differentially expressed genes between pairs of cell types
to.plot <- dt[,.(sum_genes=sum(sig,na.rm=T)),by=c("groupA","groupB")] %>%
  .[,fraction_genes:=sum_genes/total_genes] %>%
  .[,groupA:=factor(groupA,levels=opts$groups)] %>%
  .[,groupB:=factor(groupB,levels=opts$groups)]

ggplot(to.plot, aes(x=groupA, y=groupB)) +
  geom_tile(aes(fill=fraction_genes)) +
  labs(x="", y="") +
  scale_fill_gradient(low = "white", high = "red") +
  theme_classic() +
  theme(
    axis.text.y = element_text(colour="black",size=rel(1.0)),
    axis.text.x = element_text(colour="black",size=rel(1.0), angle=90, hjust=1, vjust=0.5),
    axis.ticks = element_blank()
  )

##########
## TEST ##
##########

to.plot %>% setorder(sum_genes) %>% head(n=25)

to.plot[groupB=="PGC" | groupA=="PGC"]

foo <- dt[groupA=="Mixed_mesoderm" & groupB=="PGC" & sig==T] %>% setorder(-log_padj_fdr)
foo <- dt[groupA=="Visceral_endoderm" & groupB=="Primitive_Streak" & sig==T] %>% setorder(-log_padj_fdr)

Epiblast

to.plot[groipA]