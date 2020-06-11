#####################
## Define settings ##
#####################

io <- list(); opts <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$standard.mnn.script <- "/Users/ricard/gastrulation10x/iterative_mapping/mnn/run/standard_mnn.R"
  io$iterative.mnn.script <- "/Users/ricard/gastrulation10x/iterative_mapping/mnn/run/iterative_mnn.R"
} else if(grepl("ebi",Sys.info()['nodename'])){
  io$standard.mnn.script <- "/homes/ricard/gastrulation10x/iterative_mapping/mnn/run/standard_mnn.R"
  io$iterative.mnn.script <- "/homes/ricard/gastrulation10x/iterative_mapping/mnn/run/iterative_mnn.R"
  io$tmpdir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/results/iterative_mapping/mnn/tmp"
} 

# Stages
args$stages <- c(
  "E6.5",
  # "E6.75",
  # "E7.0",
  # "E7.25",
  "E7.5",
  # "E7.75"
  # "E8.0",
  # "E8.25",
  "E8.5"
  # "mixed_gastrulation"
)

# Batches
args$test_samples <- c(
  # E7.5
  # "2",
  # "3",
  "4",
  "6"
  # "19",
  # "20"
)

# Test mode (subset cells)?
opts$test <- TRUE


#################################
## Run each plate individually ##
#################################

for (i in opts$test_samples) {
  # LSF
  if (grepl("ricard",Sys.info()['nodename'])) {
    lsf <- ""
  } else {
    lsf <- sprintf("bsub -M 30000 -n 1 -q research-rh74 -o %s/%s.txt", io$tmpdir, i)
  }

  # Run standard MNN
  cmd <- sprintf("%s Rscript %s --test_samples %s --stages %s", lsf, io$standard.mnn.script, i, paste(opts$stages, collapse=" "))
  if (isTRUE(opts$test)) cmd <- paste0(cmd, " --test")
  system(cmd)

  # Run tree-guided MNN
  cmd <- sprintf("%s Rscript %s --test_samples %s --stages %s", lsf, io$iterative.mnn.script, i, paste(opts$stages, collapse=" "))
  if (isTRUE(opts$test)) cmd <- paste0(cmd, " --test")
  system(cmd)
}


#############################
## Run all plates together ##
#############################

# LSF
if (grepl("ricard",Sys.info()['nodename'])) {
  lsf <- ""
} else {
  lsf <- sprintf("bsub -M 30000 -n 1 -q research-rh74 -o %s/all_plates.txt", io$tmpdir)
}

# Run standard MNN
cmd <- sprintf("%s Rscript %s --test_samples %s --stages %s", lsf, io$standard.mnn.script, paste(opts$test_samples, collapse=" "), paste(opts$stages, collapse=" "))
if (isTRUE(opts$test)) cmd <- paste0(cmd, " --test")
system(cmd)

# Run tree-guided MNN
cmd <- sprintf("%s Rscript %s --test_samples %s --stages %s", lsf, io$iterative.mnn.script, paste(opts$test_samples, collapse=" "), paste(opts$stages, collapse=" "))
if (isTRUE(opts$test)) cmd <- paste0(cmd, " --test")
system(cmd)
