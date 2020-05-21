# marker_genes.dt <- fread("/Users/ricard/data/gastrulation10x/results/differential/Epiblast_vs_Primitive_Streak.txt.gz") %>%
#   .[sig==T] %>%
#   .[,group:=ifelse(logFC>0,"Primitive_Streak","Epiblast")] %>%
#   .[,c("group","ens_id")]