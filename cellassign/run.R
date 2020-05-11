
################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation10x/settings.R")
  io$script <- "/Users/ricard/gastrulation10x/cellassign/cellassign.R"
} else {
  source("/homes/ricard/gastrulation10x/settings.R")
  io$script <- "/homes/ricard/gastrulation10x/cellassign/cellassign.R"
  io$tmpdir <- paste0(io$basedir,"/results/cellassign/tmp"); dir.create(io$tmpdir)
}
io$outdir <- paste0(io$basedir,"/results/cellassign")

####################
## Define options ##
####################

opts$stages <- c(
  "E6.5",
  "E6.75",
  "E7.0",
  "E7.25",
  "E7.5",
  "E7.75",
  "E8.0",
  "E8.25",
  "E8.5",
  "mixed_gastrulation"
)

# Test mode (subsetting cells)?
opts$test_mode <- TRUE

#########
## Run ##
#########

for (i in opts$stages) {
  
  # Define LSF command
  if (grepl("ricard",Sys.info()['nodename'])) {
    lsf <- ""
  } else if (grepl("ebi",Sys.info()['nodename'])) {
    lsf <- sprintf("bsub -M 90000 -n 1 -o %s/%s.txt", io$tmpdir,paste(i,collapse=" "))
  }
  cmd <- sprintf("%s Rscript %s --stages %s --outdir %s", lsf, io$script, paste(i,collapse=" "),io$outdir)
  if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test")
  
  # Run
  print(cmd)
  system(cmd)
}
