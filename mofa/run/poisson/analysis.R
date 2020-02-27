library(data.table)
library(purrr)
library(ggplot2)
library(ggpubr)

source("/Users/ricard/gastrulation10x/mofa/load_settings.R")

io <- list()
io$model.dir <- "/Users/ricard/data/gastrulation10x_mofa/hdf5/counts"
io$outdir <- "/Users/ricard/data/gastrulation10x_mofa/pdf/poisson_vs_gaussian"

##########
## Load ##
##########

model.poisson <- readRDS(paste0(io$model.dir,"/poisson_v2.rds"))
model.gaussian <- readRDS(paste0(io$model.dir,"/gaussian_v2.rds"))

###########
## Parse ##
###########

# Remove intercept factor from Poisson
model.poisson <- subset_factors(model.poisson, 2:model.poisson@dimensions[["K"]])

# Subset factors by var explained
model.gaussian <- subset_factors(model.gaussian, which(model.gaussian@cache$variance_explained$r2_per_factor[[1]][,1] > 0.01))
model.poisson <- subset_factors(model.poisson, which(model.poisson@cache$variance_explained$r2_per_factor[[1]][,1] > 0.01))

cor( colSums(model.poisson@data$view_1$group1), get_factors(model.poisson)[[1]])
cor( colSums(model.gaussian@data$view_1$group1), get_factors(model.gaussian)[[1]])

# plot_factor(model.poisson, factors = 5, color_by = "lineage_original", group_by = "lineage")
# plot_factor(model.gaussian, factors = 4, color_by = "lineage_original", group_by = "lineage")

##########
## Plot ##
##########

# plot correlation between factors
plot_factor_cor(model.poisson)
plot_factor_cor(model.gaussian)

# plot variance explained
r2.poisson.dt <- model.poisson@cache$variance_explained$r2_per_factor[[1]] %>%
  as.data.table %>% .[,factor:=as.factor(1:model.poisson@dimensions$K)] %>%
  melt(id.vars=c("factor"),variable.name="view", value.name = "r2") %>%
  .[,r2:=(r2*100)] %>% .[,cum_r2:=cumsum(r2), by="view"] %>%
  .[,model:="Poisson"]
r2.gaussian.dt <- model.gaussian@cache$variance_explained$r2_per_factor[[1]] %>%
  as.data.table %>% .[,factor:=as.factor(1:model.gaussian@dimensions$K)] %>%
  melt(id.vars=c("factor"),variable.name="view", value.name = "r2") %>%
  .[,r2:=r2*100] %>% .[,cum_r2:=cumsum(r2), by="view"] %>%
  .[,model:="Gaussian"]
r2.dt <- rbind(r2.poisson.dt,r2.gaussian.dt)

p <- ggline(r2.dt, x="factor", y="cum_r2", color="model") +
  labs(x="Factor number", y="Cumulative variance explained (%)") +
  # theme(legend.title = element_blank())
  theme(legend.title = element_blank(), legend.position = "top")

pdf(paste0(io$outdir,"/r2_vs_factor.pdf"), width=6, height=4, useDingbats = F)
print(p)
dev.off()

# compare factors
pdf(paste0(io$outdir,"/compare_factors.pdf"), width=7, height=6, useDingbats = F)
compare_factors(list("gaussian"=model.gaussian, "poisson"=model.poisson))
dev.off()

# plot dimensionality reduction
set.seed(42)
model.poisson <- run_tsne(model.poisson, perplexity = 50)
model.gaussian <- run_tsne(model.gaussian, perplexity = 50)

p.poisson <- plot_dimred(model.poisson, method="TSNE", color_by = "lineage", dot_size = 1) +
  scale_color_manual(values=opts$colors) +
  labs(title="Poisson") +
  theme(
    plot.title = element_text(hjust = 0.5, size=rel(1.2)), legend.position = "none", axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank()
  )

p.gaussian <- plot_dimred(model.gaussian, method="TSNE", color_by = "lineage", dot_size = 1) +
  scale_color_manual(values=opts$colors) +
  labs(title="Gaussian") +
  theme(
    plot.title = element_text(hjust = 0.5, size=rel(1.2)), legend.position = "none", axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank()
  )

pdf(paste0(io$outdir,"/tsne.pdf"), width=7, height=4, useDingbats = F)
cowplot::plot_grid(plotlist = list(p.gaussian,p.poisson))
dev.off()

# plot time for convergence
to.plot <- rbind(
  data.table(
    iteration = 1:length(model.poisson@training_stats$time),
    time = model.poisson@training_stats$time,
    model = "Poisson"
  ),
  data.table(
    iteration = 1:length(model.gaussian@training_stats$time),
    time = model.gaussian@training_stats$time,
    model = "Gaussian"
  )
)[time<200]

to.plot[,time:=cumsum(time), by="model"]

to.plot2 <- to.plot[,.(time=max(time)), by="model"]


p <- ggbarplot(to.plot2, x="model", y="time", fill="model") +
  labs(x="", y="Time for convergence (min)") +
  theme(
    legend.position = "none"
  )

pdf(paste0(io$outdir,"/time.pdf"), width=4, height=5, useDingbats = F)
print(p)
dev.off()

  
