
#####################
## Define settings ##
#####################

source("/Users/ricard/gastrulation10x/settings.R")
# source("/Users/ricard/gastrulation10x/iterative_mapping/mnn/stage_scores/plot_utils.R")
io$mapping.dir <- paste0(io$basedir,"/results/iterative_mapping/mnn")
io$outdir <- paste0(io$basedir,"/results/iterative_mapping/mnn/pdf")

# opts$samples <- c(
#   # E7.5
#   # "2",
#   # "3",
#   "4",
#   "6"
#   # "19",
#   # "20"
# )

opts$samples <- c("4_6")


################
## Load data  ##
################

mapping_dt <- opts$samples %>% map(function(i) {
  rbind(
    fread(sprintf("%s/%s_standard_mnn.txt.gz",io$mapping.dir,i)) %>% .[,method:="Standard MNN"],
    fread(sprintf("%s/%s_iterative_mnn.txt.gz",io$mapping.dir,i)) %>% .[,method:="Tree-guided MNN"]
    # fread(sprintf("%s/%s_standard_mnn.txt.gz",io$mapping.dir,i)) %>% .[,sample:=i] %>% .[,method:="Standard MNN"]
    # fread(sprintf("%s/%s_iterative_mnn.txt.gz",io$mapping.dir,i)) %>% .[,sample:=i] %>% .[,method:="Tree-guided MNN"]
  )
}) %>% rbindlist

mapping_dt[,.N,by="method"]
# mapping_dt[method=="Tree-guided MNN"] %>% View

############################
## Plot stage assignments ##
############################

p_list <- list()
for (i in opts$samples) {
  
  # to.plot <- mapping_dt[sample==i] %>% 
  to.plot <- mapping_dt %>% 
    merge(sample_metadata, by=c("cell")) %>% 
    .[,sample_stage:=sprintf("Sample %s (%s)", sample,stage)] %>%
    .[,.N,by=c("stage_mapped","sample_stage","method")]
  
  
  p_list[[i]] <- ggplot(to.plot, aes_string(x="stage_mapped", y="N")) +
      geom_bar(stat="identity", fill="gray70", color="black", position="dodge") +
      scale_x_discrete(drop=FALSE) + 
      facet_wrap(~sample_stage+method, scales="free_x") +
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

