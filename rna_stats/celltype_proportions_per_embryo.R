#####################
## Define settings ##
#####################

source("/Users/ricard/gastrulation10x/settings.R")
io$outfile <- paste0(io$basedir,"/results/general_stats/celltype_proportions_noExE.txt.gz")

# Remove some lineages
sample_metadata <- sample_metadata[!celltype%in%c("ExE_ectoderm","ExE_endoderm","Parietal_endoderm","Visceral_endoderm")]

dt <- sample_metadata %>% copy %>%
  .[,N:=.N,by="sample"] %>%
  .[,.(celltype_proportion=.N/unique(N)),by=c("sample","celltype")] %>%
  setorder(sample)

##########
## Save ##
##########

fwrite(dt, io$outfile, sep="\t")

