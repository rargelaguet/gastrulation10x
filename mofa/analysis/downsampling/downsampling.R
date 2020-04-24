
library(MOFA2)
library(purrr)


## Define settings ## 

io <- list()
# io$model.dir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x_mofa/hdf5/downsampling"
# io$outdir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x_mofa/pdf"
io$model.dir <- "/Users/ricard/data/gastrulation10x_mofa/hdf5/downsampling"
io$outdir <- "/Users/ricard/data/gastrulation10x_mofa/pdf"

opts <- list()
opts$downsampling_fraction <- sprintf("%.3f",seq(0, 0.25, by=0.025))

## Load precomputed models ## 

MOFAlist <- list()
for (i in opts$downsampling_fraction) {
  outfile <- sprintf("%s/model_downsample_%s.hdf5",io$model.dir,i)
  if (file.exists(outfile)) {
    MOFAlist[[i]] <- load_model(outfile, load_data = F, remove_inactive_factors = F)
    # Subset factors by variance explained
    factors <- sapply(MOFAlist[[i]]@cache$variance_explained$r2_per_factor, function(x) x[,"RNA"]>0.01)
    MOFAlist[[i]] <- subset_factors(MOFAlist[[i]], which(apply(factors,1,sum) >= 1))
  } else {
    print(sprintf("%s does not exist",outfile))
  }
}


## Assess robustness of factors ## 

# pdf(paste0(io$outdir,"/robustness.pdf"), height=4.5, width=6)
compare_factors(MOFAlist, show_rownames=F, show_colnames=F)
# dev.off()

