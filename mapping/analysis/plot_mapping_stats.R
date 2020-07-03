#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation10x/settings.R")
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/gastrulation10x/settings.R")
} 
io$mapping.dir <- paste0(io$basedir,"/results/mapping/stages")

# Cell types
opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)

# Stages
opts$stages <- c(
  "E6.5",
  "E6.75",
  "E7.0",
  "E7.25",
  "E7.5",
  "E7.75",
  "E8.0",
  "E8.25",
  "E8.5"
  # "mixed_gastrulation"
)

# Batches
opts$samples <- c(
  # they are all E7.5
  # "2",
  # "3",
  # "4",
  "6",
  "19",
  "20"
)


###################
## Load metadata ##
###################

sample_metadata <- sample_metadata %>%
  .[sample%in%opts$samples] %>%
  .[,stage:=factor(stage,levels=opts$stages)]
table(sample_metadata$sample)
table(sample_metadata$stage)

################
## Load data  ##
################

# Load mapping per sample (all cell types together)
# mapping_dt <- opts$batches %>% map(function(i) {
#   fread(sprintf("%s/mapping_mnn_%s.txt.gz",io$mapping.dir,i)) %>% .[,sample:=i]
# }) %>% rbindlist %>% .[,class:="All together"]

# Load mapping per celltype (all cell types together)
mapping_dt <- opts$samples %>% map(function(i) { 
    opts$celltypes %>% map(function(j) {
      file <- sprintf("%s/mapping_mnn_%s_%s.txt.gz",io$mapping.dir,i,j)
      if (file.exists(file)) fread(file)# %>% .[,celltype:=j]
    }) %>% rbindlist %>% .[,sample:=i]
  }) %>% rbindlist %>% .[,class:="Per cell type"] %>% .[,c("celltype_mapped","celltype_score"):=NULL]

mapping_dt %>% .[,stage_mapped:=factor(stage_mapped,levels=opts$stages)]

############################
## Plot stage assignments ##
############################

p_list <- list()
for (i in opts$samples) {
  
  to.plot <- mapping_dt[sample==i] %>%
  # to.plot <- mapping_dt %>% 
    merge(sample_metadata, by="cell") %>% 
    .[,sample_stage:=sprintf("Sample %s (%s)", i,stage)] %>%
    .[,.N,by=c("stage_mapped","sample_stage","class")]
  
  # p_list[[i]] <- ggplot(to.plot, aes_string(x="stage_mapped", y="N")) +
  ggplot(to.plot, aes_string(x="stage_mapped", y="N")) +
    geom_bar(stat="identity", fill="gray70", color="black", position="dodge") +
    scale_x_discrete(drop=FALSE) + 
    facet_wrap(~sample_stage+class, scales="free_x") +
    labs(y="Number of cells") +
    coord_flip() +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(1.3)),
      axis.title.x = element_text(color="black", size=rel(1.1)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1.3), color="black"),
      axis.text.x = element_text(size=rel(1.1), color="black")
    )
  
  # pdf(sprintf("%s/%s.pdf",io$outdir,i), width=8, height=6.5, useDingbats = F)
  print(p_list[[i]])
  # dev.off()
}

# cowplot::plot_grid(plotlist=p_list)


###############################
## Plot stage mapping scores ##
###############################

mapping_dt[,mean(stage_score),by="method"]

to.plot <- mapping_dt %>% merge(sample_metadata,by="cell") %>% 
  .[,sample_stage:=sprintf("Sample %s (%s)", sample,stage)]

ggplot(to.plot, aes_string(x="stage_score", fill="method")) +
  geom_histogram(alpha=0.5) +
  facet_wrap(~sample_stage, scales="free_y") +
  labs(x="Mapping score", y="") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(color="black", size=rel(1.2)),
    axis.text = element_text(size=rel(1), color="black")
  )





to.test <- mapping_dt %>%
  merge(sample_metadata[,c("cell","stage","celltype")]) %>%
  dcast(cell+stage+celltype~method, value.var=c("celltype_mapped","stage_mapped"))

