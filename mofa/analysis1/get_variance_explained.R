library(BioFAMtools)
library(data.table)
library(purrr)

## Define I/O ##
io <- list()
io$indir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/mofa/hdf5/stochastic"
io$outdir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/mofa/stats"

## Define options ##
opts <- list()
opts$batch_size <- c( "0.25", "0.5", "1.0" )
opts$tau <- c( "0.05", "0.15", "0.25", "0.5", "0.75", "1.0" )
opts$forgetting_rate <- c( "0.5", "0.25", "0" )
opts$ntrials <- 1


token = 1
stats <- list()

# i="0.25"; j="0.5"; k="0.75"; trial=1

# Stochastic
for (i in opts$batch_size) {
  for (j in opts$forgetting_rate) {
    for (k in opts$tau) {
      for (trial in 1:opts$ntrials) {
        print(sprintf("%s,%s,%s,%d",i,j,k,trial))
        outfile <- sprintf("%s/E6.5_k15_%s_%s_%s_%d.hdf5",io$indir,i,j,k,trial)
        
        tryCatch( {
          model <- load_model(outfile, load_training_data=TRUE, sort_factors=FALSE, set_names=TRUE)
          tmp <- calculate_variance_explained(model)
          r2 <- mean(unlist(tmp$r2_total))
          print(r2)
          stats[[token]] <- data.table(batch_size=i, forgetting_rate=j, tau=k, trial=trial, r2=r2)
        }, error = function(x) { print("Model not found") } )
        token <- token + 1
      }
    }
  }
}

# No stochastic
for (trial in 1:opts$ntrials) {
  print(sprintf("nostochastic_%d",trial))
  outfile <- sprintf("%s/nostochastic_%d.hdf5",io$indir,trial)
  
  tryCatch( {
    model <- load_model(outfile, load_training_data=TRUE, sort_factors=FALSE, set_names=TRUE)
    tmp <- calculate_variance_explained(model)
    r2 <- mean(unlist(tmp$r2_total))
    print(r2)
    stats[[token]] <- data.table(batch_size="nostochastic", forgetting_rate="nostochatic", tau="nostochastic", trial=trial, r2=r2)
  }, error = function(x) { print("Model not found") } )
  
  token <- token + 1
}

stats <- rbindlist(stats)
fwrite(stats, paste0(io$outdir,"/r2_nostochastic.txt"), col.names=T, quote=F, sep="\t")