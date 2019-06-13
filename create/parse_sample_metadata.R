library(data.table)
library(purrr)

sample_metadata <- fread("/Users/ricard/data/gastrulation10x/original/meta.tab") %>%
  .[,lineage10x:=celltype] %>%
  .[,lineage10x_2:=lineage10x]

foo <- sample_metadata[stage%in%c("E6.5","E6.75","E7.0","E7.25","E7.5","mixed_gastrulation")] %>%
  # Mesoderm
  .[lineage10x%in%c("Caudal Mesoderm","Somitic mesoderm","Pharyngeal mesoderm","Paraxial mesoderm","ExE mesoderm","Mesenchyme","Intermediate mesoderm", "Mixed mesoderm", "Nascent mesoderm","Allantois"), lineage10x_2:="Mesoderm"] %>%
  .[lineage10x%in%c("Blood progenitors 2", "Blood progenitors 1", "Endothelium", "Erythroid2", "Erythroid3", "Haematoendothelial progenitors","Erythroid1","Cardiomyocytes"), lineage10x_2:="Blood"] %>%
  # Endoderm
  .[lineage10x%in%c("Gut","Def. endoderm","Notochord"), lineage10x_2:="Endoderm"] %>%
  .[lineage10x%in%c("Parietal endoderm","Visceral endoderm"), lineage10x_2:="ExE endoderm"] %>%
  # Primitive streak
  .[lineage10x%in%c("Caudal epiblast","Anterior Primitive Streak"), lineage10x_2:="Primitive Streak"] %>%
  # Ectoderm
  .[lineage10x%in%c("NMP","Spinal cord", "Forebrain/Midbrain/Hindbrain", "Caudal neurectoderm","Rostral neurectoderm","Surface ectoderm"), lineage10x_2:="Ectoderm"]
  # Brain

bar <- sample_metadata[!stage%in%c("E6.5","E6.75","E7.0","E7.25","E7.5","mixed_gastrulation")]

sample_metadata <- rbind(foo,bar)

fwrite(sample_metadata, "/Users/ricard/data/gastrulation10x/sample_metadata2.txt", sep="\t", col.names=T, row.names=F, na="NA", quote=F)



###

table(sample_metadata[stage=="E6.5",lineage10x_2])
