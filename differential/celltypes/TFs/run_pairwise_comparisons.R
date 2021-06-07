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
io$outdir <- paste0(io$basedir,"/results/differential/TFs"); dir.create(io$outdir, showWarnings = F)
io$tmpdir <- paste0(io$basedir,"/tmp"); dir.create(io$outdir, showWarnings = F)

#############
## Options ##
#############

# Testing mode
opts$test_mode <- FALSE

# Define cell types
# opts$groups <- opts$celltypes
opts$groups <- names(which(table(sample_metadata[stage%in%opts$stages,celltype])>=100))

#########
## Run ##
#########

for (i in 1:length(opts$groups)) {
  groupA <- opts$groups[[i]]
  for (j in 1:length(opts$groups)) {
    if (i!=j) {
      groupB <- opts$groups[[j]]
      outfile <- sprintf("%s/%s_vs_%s.txt.gz", io$outdir,groupA,groupB)
      
      # Define LSF command
      if (grepl("ricard",Sys.info()['nodename'])) {
        lsf <- ""
      } else if (grepl("ebi",Sys.info()['nodename'])) {
        lsf <- sprintf("bsub -M 18000 -n 1 -o %s/%s_vs_%s.txt", io$tmpdir,groupA,groupB)
      }
      cmd <- sprintf("%s Rscript %s --stages %s --groupA %s --groupB %s --test %s --outfile %s", lsf, io$script, paste(opts$stages, collapse=" "), groupA, groupB, test, outfile)
      if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test_mode")
      
      # Run
      print(cmd)
      system(cmd)
    }
  }
}

