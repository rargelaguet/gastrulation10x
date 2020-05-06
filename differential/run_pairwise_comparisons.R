#########
## I/O ##
#########

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation10x/settings.R")
  io$script <- "/Users/ricard/gastrulation10x/differential/differential.R"
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/gastrulation10x/settings.R")
  io$script <- "/homes/ricard/gastrulation10x/differential/differential.R"
  io$tmpdir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/results/differential/tmp"; dir.create(io$tmpdir, showWarnings=F)
} else {
  stop("Computer not recognised")
}
io$outdir <- paste0(io$basedir,"/results/differential"); dir.create(io$outdir, showWarnings = F)

#############
## Options ##
#############

# Define cell types
opts$groups <- opts$celltypes.1

# Statistical test
# opts$test <- c("edgeR","t-test")
opts$statistical.test <- c("edgeR")

# Testing mode
opts$test_mode <- FALSE

for (test in opts$statistical.test) {
  for (i in 1:(length(opts$groups)-1)) {
    groupA <- opts$groups[[i]]
    for (j in 1:length(opts$groups)) {
      if (i!=j) {
        groupB <- opts$groups[[j]]
        outfile <- sprintf("%s/%s_vs_%s.txt.gz", io$outdir,groupA,groupB)# %>% 
          # stringr::str_replace_all(.," ","-")
        
        # Define LSF command
        if (grepl("ricard",Sys.info()['nodename'])) {
          lsf <- ""
        } else if (grepl("ebi",Sys.info()['nodename'])) {
          lsf <- sprintf("bsub -M 15000 -n 1 -q research-rh74 -o %s/%s_vs_%s.txt", io$tmpdir,groupA,groupB)
        }
        cmd <- sprintf("%s Rscript %s --groupA %s --groupB %s --test %s --outfile %s", lsf, io$script, groupA, groupB, test, outfile)
        if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test_mode")
        
        # Run
        print(cmd)
        system(cmd)
      }
    }
  }
}
