library(data.table)
library(purrr)
library(ggpubr)

embryonic.endoderm <- c("Def. endoderm","Gut", "Parietal endoderm")
ExE.endoderm <- c("ExE endoderm","Visceral endoderm", "Parietal endoderm")
ExE.ectoderm <- c("ExE ectoderm")
ExE.mesoderm <- c("ExE mesoderm")


sample_metadata <- fread("/Users/ricard/data/gastrulation10x/sample_metadata2.txt") %>%
  .[,sample:=as.factor(sample)]
sample_metadata[,ExE:=lineage10x%in%ExE.endoderm]
sample_metadata[,embryonic:=lineage10x%in%embryonic.endoderm]

to.plot <- sample_metadata[,.(ratio=sum(ExE)/sum(embryonic)), by=c("sample","stage")] %>%
  .[stage%in%c("E7.25","E7.5","E7.75","E8.0","E8.25","E8.5")]

ggbarplot(to.plot, x="sample", y="ratio", fill="gray70") +
  facet_wrap(~stage, scales = "free_x", nrow=2) +
  theme(
    axis.text.x = element_blank()
  )
