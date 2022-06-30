#########
## I/O ##
#########

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation10x/settings.R")
  io$script <- "/Users/ricard/gastrulation10x/differential/celltypes/TFs/differential.R"
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/gastrulation10x/settings.R")
  io$script <- "/homes/ricard/gastrulation10x/differential/celltypes/TFs/differential.R"
} else {
  stop("Computer not recognised")
}
io$outdir <- paste0(io$basedir,"/results/differential/celltypes/TFs"); dir.create(io$outdir, showWarnings = F)
io$tmpdir <- paste0(io$outdir,"/tmp"); dir.create(io$tmpdir, showWarnings = F)

#############
## Options ##
#############

# Testing mode
opts$test_mode <- FALSE

# Define cell types
opts$groups <- opts$celltypes# %>% head(n=3)

#########
## Run ##
#########

for (i in 1:length(opts$groups)) {
  groupA <- opts$groups[[i]]
  for (j in i:length(opts$groups)) {
    if (i!=j) {
      groupB <- opts$groups[[j]]
      outfile <- sprintf("%s/%s_vs_%s.txt.gz", io$outdir,groupA,groupB)
      
      # Define LSF command
      if (grepl("ricard",Sys.info()['nodename'])) {
        lsf <- ""
      } else if (grepl("ebi",Sys.info()['nodename'])) {
        lsf <- sprintf("bsub -M 16000 -n 1 -o %s/%s_vs_%s.txt", io$tmpdir,groupA,groupB)
      }
      cmd <- sprintf("%s Rscript %s --groupA %s --groupB %s --outfile %s", lsf, io$script, groupA, groupB, outfile)
      if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test_mode")
      
      # Run
      print(cmd)
      system(cmd)
    }
  }
}

