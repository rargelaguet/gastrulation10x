#########
## I/O ##
#########

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation10x/settings.R")
  io$script <- "/Users/ricard/gastrulation10x/differential/differential.R"
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/gastrulation10x/settings.R")
  io$script <- "/homes/ricard/gastrulation10x/differential/rna/differential.R"
  io$tmpdir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/results/differential/tmp"; dir.create(io$tmpdir, showWarnings=F)
} else {
  stop("Computer not recognised")
}
io$outdir <- paste0(io$basedir,"/results/differential"); dir.create(io$outdir, showWarnings = F)

#############
## Options ##
#############

opts <- list()

# Define cell types
opts$groups <- c(
  "Epiblast",
  "Primitive Streak",
  "ExE ectoderm",
  "Visceral endoderm"
  # "ExE endoderm",
  # "Nascent mesoderm",
  # "Neurectoderm",
  # "Blood progenitors",
  # "Mixed mesoderm",
  # "ExE mesoderm",
  # "Pharyngeal mesoderm",
  # "Caudal epiblast",
  # "PGC",
  # "Mesenchyme",
  # "Haematoendothelial progenitors",
  # "Surface ectoderm",
  # "Gut",
  # "Paraxial mesoderm",
  # "Notochord",
  # "Somitic mesoderm",
  # "Caudal Mesoderm",
  # "Erythroid",
  # "Def. endoderm",
  # "Parietal endoderm",
  # "Allantois",
  # "Anterior Primitive Streak",
  # "Endothelium",
  # "Forebrain_Midbrain_Hindbrain",
  # "Spinal cord",
  # "Cardiomyocytes",
  # "NMP",
  # "Neural crest"
) %>% stringr::str_replace_all(.," ","-")

# Statistical test
# opts$test <- c("edgeR","t-test")
opts$test <- c("edgeR")

# Testing mode
opts$test_mode <- TRUE

for (test in opts$test) {
  for (i in 1:(length(opts$groups)-1)) {
    groupA <- opts$groups[[i]]
    for (j in (i+1):length(opts$groups)) {
      groupB <- opts$groups[[j]]
      outfile <- sprintf("%s/%s_vs-%s.txt.gz", io$outdir,groupA,groupB) %>% 
        stringr::str_replace_all(.," ","-")
      
      # Define LSF command
      if (grepl("ricard",Sys.info()['nodename'])) {
        lsf <- ""
      } else if (grepl("ebi",Sys.info()['nodename'])) {
        lsf <- sprintf("bsub -M 3000 -n 1 -q research-rh74 -o %s/%s_vs_%s.txt", io$tmpdir,groupA,groupB)
      }
      cmd <- sprintf("%s Rscript %s --groupA %s --groupB %s --test %s --outfile %s", lsf, io$script, groupA, groupB, test, outfile)
      if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test_mode")
      
      # Run
      print(cmd)
      system(cmd)
    }
  }
}
