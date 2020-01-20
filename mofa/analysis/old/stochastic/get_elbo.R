# library(BioFAMtools)
# devtools::load_all("/homes/ricard/biofam/BioFAMtools")
library(data.table)
library(purrr)
library(rhdf5)

## Define I/O ##
io <- list()
io$indir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/mofa/hdf5"
io$outdir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/mofa/stats"

## Define options ##
opts <- list()
opts$batch_size <- c( "0.25", "0.5", "1.0" )
opts$tau <- c( "0.05", "0.15", "0.25", "0.5", "0.75", "1.0" )
opts$forgetting_rate <- c( "0.5", "0.25", "0" )
opts$ntrials <- 1


token = 1
stats <- list()

# Stochastic
for (i in opts$batch_size) {
  for (j in opts$forgetting_rate) {
    for (k in opts$tau) {
      for (trial in 1:opts$ntrials) {
        print(sprintf("%s,%s,%s,%d",i,j,k,trial))
        outfile <- sprintf("%s/E6.5_k15_%s_%s_%s_%d.hdf5",io$indir,i,j,k,trial)
        
        tryCatch( {
          # model <- load_model(outfile, load_training_data=FALSE, sort_factors=FALSE, set_names=FALSE)
          training_stats <- h5read(outfile, 'training_stats', read.attributes=T)
          elbo = training_stats$elbo
          time = training_stats$time
          stats[[token]] <- data.table(iter=1:length(elbo), batch_size=i, forgetting_rate=j, tau=k, trial=trial, elbo=elbo, time=cumsum(time))
        }, error = function(x) { print("Model not found") } )
        token <- token + 1
    }
    }
  }
}

# No stochastic
for (k in 1:opts$ntrials) {
  print(sprintf("No stochastic %d",k))
  outfile <- sprintf("%s/nostochastic_%d.hdf5",io$indir,k)
  
  tryCatch( {
    # model <- load_model(outfile, load_training_data=FALSE, sort_factors=FALSE, set_names=FALSE)
    training_stats <- h5read(outfile, 'training_stats', read.attributes=T)
    elbo = training_stats$elbo
    time = training_stats$time
    stats[[token]] <- data.table(iter=1:length(elbo), tau="nostochastic", forgetting_rate="nostochastic", batch_size="nostochastic", trial=k, elbo=elbo, time=cumsum(time))
  }, error = function(x) { print("Model not found") } )
  
  token <- token + 1
}


stats <- rbindlist(stats)

fwrite(stats, paste0(io$outdir,"/stats.txt"), col.names=T, quote=F, sep="\t")