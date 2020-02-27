library(MOFA2)
library(data.table)
library(purrr)
library(ggpubr)

#####################
## Define settings ##
#####################

io <- list()
io$model <- "/Users/ricard/data/gastrulation10x_mofa/hdf5/model_1.hdf5"
io$outdir <- "/Users/ricard/data/gastrulation10x_mofa/pdf/variance_explained"

################
## Load model ##
################

model <- load_model(io$model)

factors(model) <- as.character(1:10)

plot_factor_cor(model)

#######################################
## Plot variance explained vs factor ##
#######################################

to.plot <- get_variance_explained(model, as.data.frame = T)[["r2_per_factor"]] %>%
  as.data.table %>% setnames("value","r2") %>%
  .[,r2:=r2*100] %>%
  .[,cum_r2:=cumsum(r2), by="group"]

p <- ggline(to.plot, x="factor", y="cum_r2", color="group") +
  labs(x="Factor number", y="Variance explained (%)") +
  # coord_cartesian(ylim=c(10,55)) +
  scale_color_brewer(palette="Dark2") +
  theme(legend.title = element_blank())

pdf(paste0(io$outdir,"/r2_vs_factor.pdf"), width=6, height=4, useDingbats = F)
print(p)
dev.off()
