
library(MOFA2)
library(purrr)


## Define settings ## 

io <- list()
io$model.dir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x_mofa/hdf5/downsampling"
io$outdir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x_mofa/pdf"

opts <- list()
opts$downsampling_fraction <- sprintf("%.3f",seq(0, 0.25, by=0.025))


## Load precomputed models ## 

MOFAlist <- list()
for (i in opts$downsampling_fraction) {
  outfile <- sprintf("%s/model_downsample_%s.hdf5",io$model.dir,i)
  MOFAlist[[i]] <- load_model(outfile)
}


## Assess robustness of factors ## 

pdf(paste0(io$outdir,"/robustness.pdf"), height=7, width=8.5)
compare_factors(MOFAlist, show_rownames=F, show_colnames=F)
dev.off()

