source("/Users/ricard/gastrulation10x/settings.R")

library(ggrepel)

#####################
## I/O and options ##
#####################

io$outdir <- paste0(io$basedir,"/results/blood_trajectory")

names(opts$celltype.colors) <- names(opts$celltype.colors) %>% stringr::str_replace_all("_", " ")

opts$haem_colours = c(
  "Mes1"= "#c4a6b2",
  "Mes2"= "#ca728c",
  
  "Cardiomyocytes" =  "#B51D8D",  
  
  "BP1" = "#6460c5",
  "BP2" = "#96b8e4",
  "Haem3"= "#02f9ff",
  "BP3" = "#07499f",
  "BP4" = "#036ef9",
  
  "Haem1"= "#bb22a7",
  "Haem2" = "#f695e9",
  "Haem4" = "#4c4a81",
  
  "EC1"= "#006737",
  
  "EC2" = "#5eba46",
  "EC3" = "#818068",
  "EC4"="#d6de22",
  "EC5"="#5f7238",
  "EC6"="#2ab572",
  "EC7"="#000000",
  "EC8"="#a0cb3b",
  
  "Ery1"="#f67a58",
  "Ery2" ="#a26724",
  "Ery3"="#cdaf7f",
  "Ery4"= "#625218",
  
  "My" = "#c62127",
  "Mk"= "#f6931d"
)

###############
## Load data ##
###############

# Load atlas metadata (which includes precomputed UMAP coordinates)
meta_atlas <- fread("/Users/ricard/data/gastrulation10x/original/sample_metadata.txt.gz") %>%
  .[!is.na(haem_gephiY)] %>%
  .[,c("cell","sample", "stage", "celltype", "umapX", "umapY", "haem_gephiX", "haem_gephiY","haem_subclust")] %>%
  droplevels()

################
## Parse data ##
################

# names(opts$celltype.colors) <- names(opts$celltype.colors) %>% stringr::str_replace_all("_"," ")
# opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% unique(to.plot$aggregated_celltype)]

to.plot <- meta_atlas

###################################
## Plot dimensionality reduction ##
###################################

p <- ggplot(to.plot, aes(x=haem_gephiX, y=haem_gephiY)) +
  ggrastr::geom_point_rast(aes(colour=haem_subclust), size=0.05) + scale_color_manual(values=opts$haem_colours) +
  # ggrastr::geom_point_rast(aes(colour=celltype), size=0.05) + scale_color_manual(values=opts$celltype.colors) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )

# Label clusters
groups <- unique(to.plot$haem_subclust)
labels.loc <- lapply(groups, function(i) {
    data.table(t(apply(X = to.plot[haem_subclust==i, c("haem_gephiX","haem_gephiY")], MARGIN = 2, FUN = median, na.rm = TRUE))) %>%
      .[,haem_subclust:=i] %>% return
  }) %>% rbindlist %>% setnames(c("haem_gephiX","haem_gephiY","haem_subclust"))
p <- p + geom_text_repel(aes_string(label="haem_subclust"), data=labels.loc)

pdf(paste0(io$outdir,"/blood_trajectory_haem_subclust.pdf"), width=7, height=6, useDingbats = F)
print(p)
dev.off()
