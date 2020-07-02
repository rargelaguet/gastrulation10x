#########
## I/O ##
#########

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation10x/settings.R")
  io$script <- "/Users/ricard/gastrulation10x/differential/stages/differential.R"
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/gastrulation10x/settings.R")
  io$script <- "/homes/ricard/gastrulation10x/differential/stages/differential.R"
  io$tmpdir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/results/differential/E7.5/tmp"; dir.create(io$tmpdir, showWarnings=F)
} else {
  stop("Computer not recognised")
}
io$outdir <- paste0(io$basedir,"/results/differential/stages"); dir.create(io$outdir, showWarnings = F)

#############
## Options ##
#############

# Define stages
opts$stages <- c(
  "E6.5",
  # "E6.75",
  # "E7.0"
  # "E7.25",
  "E7.5"
  # "E7.75"
  # "E8.0",
  # "E8.25",
  # "E8.5"
  # "mixed_gastrulation"
)

opts$celltypes = c(
  # "Epiblast",
  "Primitive_Streak"
  # "Caudal_epiblast",
  # "PGC",
  # "Anterior_Primitive_Streak",
  # "Notochord",
  # "Def._endoderm",
  # "Gut",
  # "Nascent_mesoderm",
  # "Mixed_mesoderm",
  # "Intermediate_mesoderm",
  # "Caudal_Mesoderm",
  # "Paraxial_mesoderm",
  # "Somitic_mesoderm",
  # "Pharyngeal_mesoderm",
  # "Cardiomyocytes",
  # "Allantois",
  # "ExE_mesoderm",
  # "Mesenchyme",
  # "Haematoendothelial_progenitors",
  # "Endothelium",
  # "Blood_progenitors_1",
  # "Blood_progenitors_2",
  # "Erythroid1",
  # "Erythroid2",
  # "Erythroid3",
  # "NMP",
  # "Rostral_neurectoderm",
  # "Caudal_neurectoderm",
  # "Neural_crest",
  # "Forebrain_Midbrain_Hindbrain",
  # "Spinal_cord",
  # "Surface_ectoderm",
  # "Visceral_endoderm",
  # "ExE_endoderm",
  # "ExE_ectoderm",
  # "Parietal_endoderm"
)

# Testing mode
opts$test_mode <- FALSE

# Update sample metadata
sample_metadata <- sample_metadata[celltype%in%opts$celltypes & stage%in%opts$stages]

#########
## Run ##
#########

for (i in opts$celltypes) {
  opts$stages <- names(which(table(sample_metadata[celltype==i,stage])>=50))
  for (j in 1:length(opts$stages)) {
    stageA <- opts$stages[[j]]
    for (k in 1:length(opts$stages)) {
      if (j!=k) {
        stageB <- opts$stages[[k]]
        outfile <- sprintf("%s/%s_%s_vs_%s.txt.gz", io$outdir,i,stageA,stageB)
        
        # Define LSF command
        if (grepl("ricard",Sys.info()['nodename'])) {
          lsf <- ""
        } else if (grepl("ebi",Sys.info()['nodename'])) {
          lsf <- sprintf("bsub -M 18000 -n 1 -o %s/%s_%s_vs_%s.txt", io$tmpdir,i,stageA,stageB)
        }
        cmd <- sprintf("%s Rscript %s --celltype %s --stageA %s --stageB %s --test edgeR --outfile %s", lsf, io$script, i, stageA, stageB, outfile)
        if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test_mode")
        
        # Run
        print(cmd)
        system(cmd)
      }
    }
  }
}
