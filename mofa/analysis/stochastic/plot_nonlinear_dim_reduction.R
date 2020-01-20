# library(BioFAMtools)
devtools::load_all("/homes/ricard/biofam/BioFAMtools")
library(umap)
library(Rtsne)
library(data.table)
library(purrr)
library(ggplot2)


## Define I/O ##
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation10x/mofa/load_settings.R")
} else {
  source("/homes/ricard/gastrulation10x/mofa/load_settings.R")
}
io$indir <- paste0(io$basedir,"/mofa/hdf5/stochastic")

## Define options ##
opts$batch_size <- c( "0.25", "0.5" )
opts$tau <- c( "0.05", "0.15", "0.25", "0.5", "0.75", "1.0" )
opts$forgetting_rate <- c( "0.5", "0.25", "0" )
opts$ntrials <- 1
opts$algorithms <- c("tsne","umap")


## Load data ##
sample_metadata <- fread(io$sample.metadata, showProgress=F)

# i=opts$batch_size[1]; j=opts$forgetting_rate[1]; k=opts$tau[1]; trial=1; algorithm="tsne"

for (i in opts$batch_size) {
  for (j in opts$forgetting_rate) {
    for (k in opts$tau) {
      for (trial in 1:opts$ntrials) {
        for (algorithm in opts$algorithms) {

          # load model
          model_file <- sprintf("%s/E6.5_k15_%s_%s_%s_%d.hdf5",io$indir,i,j,k,trial)
          print(model_file)
          model <- load_model(file=model_file, load_training_data=FALSE, sort_factors=FALSE, set_names=TRUE)

          # Parse sample metadata
          sample_metadata_filt <- sample_metadata %>% setkey(cell) %>% .[unname(unlist(samples_names(model)))]

          # subset factors
          r2 <- calculate_variance_explained(model)$r2_per_factor
          factors <- sapply(r2, function(x) x[,"RNA"]>0.001)
          
          if (length(factors)!=0) {
            model <- subset_factors(model, which(apply(factors,1,sum) >= 1))
            print(get_dimensions(model)[["K"]])
            factors_names(model) <- paste("Factor",1:get_dimensions(model)[["K"]], sep=" ")

            # fetch factors
            Z <- get_factors(model) %>% do.call("rbind",.)
            
            # do dimensionality reudction
            set.seed(1)
            if (algorithm=="tsne") {
              tsne <- Rtsne(Z, check_duplicates=FALSE, pca=FALSE, theta=0.5, dims=2)
              Z.out <- tsne$Y
            } else if (algorithm=="umap") {
              umap.defaults$n_neighbors <- 25
              umap.defaults$min_dist <- 0.5
              umap.out <- umap(Z, config = umap.defaults)
              Z.out <- umap.out$layout
            }
            
            # plot
            to.plot <- Z.out %>% as.data.table %>% .[,cell:=rownames(Z)] %>%
                merge(sample_metadata_filt, by="cell")

            p <- ggplot(to.plot, aes(x=V1, y=V2, color=`celltype`)) +
              geom_point(alpha=0.7, size=1.5) +
              labs(x="Dimension 1", y="Dimension 2") +
              scale_color_manual(values=opts$colors) +
              guides(colour = guide_legend(override.aes = list(size=3))) +
              theme(
                legend.position = "none",
                legend.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank()
              )
            
            pdf(sprintf("%s/%s_MOFA_celltype_%s_%s_%s_%d.pdf",io$pdfdir,algorithm,i,j,k,trial), width=5, height=4.5, useDingbats = F)
            print(p)
            dev.off()

          }
        }
      }
    }
  }
}
