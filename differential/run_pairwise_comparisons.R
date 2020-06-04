#########
## I/O ##
#########

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation10x/settings.R")
  io$script <- "/Users/ricard/gastrulation10x/differential/differential.R"
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/gastrulation10x/settings.R")
  io$script <- "/homes/ricard/gastrulation10x/differential/differential.R"
  io$tmpdir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/results/differential/E8.5/tmp"; dir.create(io$tmpdir, showWarnings=F)
} else {
  stop("Computer not recognised")
}
io$outdir <- paste0(io$basedir,"/results/differential"); dir.create(io$outdir, showWarnings = F)

#############
## Options ##
#############

# Define stages
opts$stages <- c(
  # "E6.5",
  # "E6.75",
  # "E7.0",
  # "E7.25",
  # "E7.5",
  # "E7.75"
  # "E8.0",
  # "E8.25",
  "E8.5"
  # "mixed_gastrulation"
)

# Statistical test
opts$statistical.test <- c("edgeR")

# Testing mode
opts$test_mode <- FALSE

# Define cell types
# opts$groups <- opts$celltypes.1
opts$groups <- names(which(table(sample_metadata[stage%in%opts$stages,celltype])>=100))

#########
## Run ##
#########

for (test in opts$statistical.test) {
  for (i in 1:length(opts$groups)) {
  # for (i in length(opts$groups)) {
    groupA <- opts$groups[[i]]
    for (j in 1:length(opts$groups)) {
      if (i!=j) {
        groupB <- opts$groups[[j]]
        outfile <- sprintf("%s/E8.5/%s_vs_%s.txt.gz", io$outdir,groupA,groupB)
        
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
}


##############################
## Run selected comparisons ##
##############################

# opts$comparisons <- list(
#   c("groupA"="Erythroid1",          "groupB"="Mixed_mesoderm"),
#   c("groupA"="Blood_progenitors_2", "groupB"="Caudal_Mesoderm"),
#   c("groupA"="Allantois",           "groupB"="Haematoendothelial_progenitors"),
#   c("groupA"="Blood_progenitors_2", "groupB"="NMP"),
#   c("groupA"="Erythroid1",          "groupB"="Anterior_Primitive_Streak"),
#   c("groupA"="Erythroid1",          "groupB"="Haematoendothelial_progenitors")
# )

# for (comparison in opts$comparisons) {
#   groupA <- comparison[["groupA"]]; groupB <- comparison[["groupB"]]
#   for (test in opts$statistical.test) {
#     outfile <- sprintf("%s/%s_vs_%s.txt.gz", io$outdir,groupA,groupB)

#     # Define LSF command
#     if (grepl("ricard",Sys.info()['nodename'])) {
#       lsf <- ""
#     } else if (grepl("ebi",Sys.info()['nodename'])) {
#       lsf <- sprintf("bsub -M 15000 -n 1 -q research-rh74 -o %s/%s_vs_%s.txt", io$tmpdir,groupA,groupB)
#     }
#     cmd <- sprintf("%s Rscript %s --groupA %s --groupB %s --test %s --outfile %s", lsf, io$script, groupA, groupB, test, outfile)
#     if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test_mode")

#     # Run
#     print(cmd)
#     system(cmd)
#   }
# }

