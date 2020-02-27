library(MOFA2)

io <- list()
io$metadata <- "/Users/ricard/data/gastrulation10x/sample_metadata.txt.gz"
io$outdir <- "/Users/ricard/data/gastrulation10x_mofa/hdf5/counts"

##################
## RUN POISSON ##
##################

# Cap values
m_counts[m_counts>=100] <- 100

# MOFAobject <- create_mofa(list(t(m_counts[[1]])))
MOFAobject <- create_mofa(list(m_counts))

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 9
model_opts$likelihoods <- "poisson"
# model_opts$ard_factors <- FALSE
# model_opts$ard_weights <- FALSE
# model_opts$spikeslab_factors <- FALSE
# model_opts$spikeslab_weights <- FALSE

train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "fast"

MOFAobject <- prepare_mofa(MOFAobject,
  model_options = model_opts,
  training_options = train_opts
)

model.poisson <- run_mofa(MOFAobject)

##################
## RUN GAUSSIAN ##
##################

MOFAobject <- create_mofa(list(m_logcounts))

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 8

train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "fast"

MOFAobject <- prepare_mofa(MOFAobject,
  model_options = model_opts,
  training_options = train_opts
)

model.gaussian <- run_mofa(MOFAobject)

#########################
## Add sample metadata ##
#########################

sample_metadata <- fread(io$metadata) %>%
  .[,c("cell","lineage10x","lineage10x_2")] %>% 
  setnames(c("sample","lineage_original","lineage")) %>%  
  .[,group:="group1"]
cells <- as.character(unname(unlist(MOFA2::samples(model.poisson))))
sample_metadata_filt <- sample_metadata[sample%in%cells] %>% setkey(sample) %>% .[cells]
stopifnot(all(cells==sample_metadata_filt$sample))
samples_metadata(model.gaussian) <- sample_metadata_filt
samples_metadata(model.poisson) <- sample_metadata_filt

##########
## Save ##
##########

saveRDS(model.poisson, paste0(io$outdir,"/poisson_v2.rds"))
saveRDS(model.gaussian, paste0(io$outdir,"/gaussian_v2.rds"))





##################
## TEST MOFA V1 ##
##################

# library(MOFA)
# 
# MOFAobject <- createMOFAobject(list(m_counts))
# 
# model_opts <- getDefaultModelOptions(MOFAobject)
# model_opts$numFactors <- 5
# 
# train_opts <- getDefaultTrainOptions()
# train_opts$tolerance <- 1.0
# 
# MOFAobject <- prepareMOFA(MOFAobject, ModelOptions = model_opts, TrainOptions = train_opts)
# 
# runMOFA(MOFAobject)
