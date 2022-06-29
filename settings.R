suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(SingleCellExperiment))

#########
## I/O ##
#########

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/gastrulation10x"
  io$gene_metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt"
} else if (grepl("rargelaguet.local",Sys.info()['nodename'])) {
  io$basedir <- "/Users/rargelaguet/data/gastrulation10x"
  io$gene_metadata <- "/Users/rargelaguet/data/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt"
} else if (grepl("ebi",Sys.info()['nodename'])) {
  io$basedir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x"
  io$gene_metadata <- "/hps/nobackup2/research/stegle/users/ricard/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt"
} else {
  stop("Computer not recognised")
}

io$metadata <- paste0(io$basedir,"/sample_metadata.txt.gz")
# io$marker_genes.stringent <- paste0(io$basedir,"/results/marker_genes/marker_genes_stringent.tsv.gz")
io$marker_genes <- paste0(io$basedir,"/results/marker_genes/all_stages/marker_genes.txt.gz")
io$average_expression_per_celltype <- paste0(io$basedir,"/results/marker_genes/all_stages/avg_expr_per_celltype_and_gene.txt.gz")
io$sce <- paste0(io$basedir,"/processed/SingleCellExperiment.rds")
io$cellassign.dir <- paste0(io$basedir,"/results/cellassign")

io$paga.connectivity <- paste0(io$basedir,"/results/paga/paga_connectivity.csv")
io$paga.coordinates <- paste0(io$basedir,"/results/paga/paga_coordinates.csv")

#############
## Options ##
#############

opts <- list()

opts$celltype.colors = c(
	"Epiblast" = "#635547",
	"Primitive_Streak" = "#DABE99",
	"Caudal_epiblast" = "#9e6762",
	"PGC" = "#FACB12",
	"Anterior_Primitive_Streak" = "#c19f70",
	"Notochord" = "#0F4A9C",
	"Def._endoderm" = "#F397C0",
	"Gut" = "#EF5A9D",
	"Nascent_mesoderm" = "#C594BF",
	"Mixed_mesoderm" = "#DFCDE4",
	"Intermediate_mesoderm" = "#139992",
	"Caudal_Mesoderm" = "#3F84AA",
	"Paraxial_mesoderm" = "#8DB5CE",
	"Somitic_mesoderm" = "#005579",
	"Pharyngeal_mesoderm" = "#C9EBFB",
	"Cardiomyocytes" = "#B51D8D",
	"Allantois" = "#532C8A",
	"ExE_mesoderm" = "#8870ad",
	"Mesenchyme" = "#cc7818",
	"Haematoendothelial_progenitors" = "#FBBE92",
	"Endothelium" = "#ff891c",
	"Blood_progenitors_1" = "#f9decf",
	"Blood_progenitors_2" = "#c9a997",
	"Erythroid1" = "#C72228",
	"Erythroid2" = "#f79083",
	"Erythroid3" = "#EF4E22",
	"NMP" = "#8EC792",
	"Rostral_neurectoderm" = "#65A83E",
	"Caudal_neurectoderm" = "#354E23",
	"Neural_crest" = "#C3C388",
	"Forebrain_Midbrain_Hindbrain" = "#647a4f",
	"Spinal_cord" = "#CDE088",
	"Surface_ectoderm" = "#f7f79e",
	"Visceral_endoderm" = "#F6BFCB",
	"ExE_endoderm" = "#7F6874",
	"ExE_ectoderm" = "#989898",
	"Parietal_endoderm" = "#1A1A1A",
	
	# Additional
	"Erythroid" = "#EF4E22",
	"Blood_progenitors" = "#c9a997",
	"Neurectoderm" = "#65A83E"
)

opts$stage.colors = c(
	"E8.5" = "#440154FF",
	"E8.25" = "#472D7BFF",
	"E8.0" = "#3B528BFF",
	"E7.75" = "#2C728EFF",
	"E7.5" = "#21908CFF",
	"E7.25" = "#27AD81FF",
	"E7.0" = "#5DC863FF",
	"E6.75" = "#AADC32FF",
	"E6.5" = "#FDE725FF"
)
# opts$celltype.colors.2 = c(
# 	"Epiblast" = "#635547",
# 	"Primitive_Streak" = "#DABE99",
# 	"Caudal_epiblast" = "#9e6762",
# 	"PGC" = "#FACB12",
# 	"Anterior_Primitive_Streak" = "#c19f70",
# 	"Notochord" = "#0F4A9C",
# 	"Def._endoderm" = "#F397C0",
# 	"Gut" = "#EF5A9D",
# 	"Nascent_mesoderm" = "#C594BF",
# 	"Mixed_mesoderm" = "#DFCDE4",
# 	# "Intermediate mesoderm" = "#139992",
# 	"Caudal_Mesoderm" = "#3F84AA",
# 	"Paraxial_mesoderm" = "#8DB5CE",
# 	"Somitic_mesoderm" = "#005579",
# 	"Pharyngeal_mesoderm" = "#C9EBFB",
# 	"Cardiomyocytes" = "#B51D8D",
# 	"Allantois" = "#532C8A",
# 	"ExE_mesoderm" = "#8870ad",
# 	"Mesenchyme" = "#cc7818",
# 	"Haematoendothelial_progenitors" = "#FBBE92",
# 	"Endothelium" = "#ff891c",
# 	# "Blood_progenitors_1" = "#f9decf",
# 	# "Blood_progenitors_2" = "#c9a997",
# 	"Blood_progenitors" = "#c9a997",
# 	# "Erythroid1" = "#C72228",
# 	# "Erythroid2" = "#f79083",
# 	# "Erythroid3" = "#EF4E22",
# 	"Erythroid" = "#EF4E22",
# 	"NMP" = "#8EC792",
# 	# "Rostral_neurectoderm" = "#65A83E",
# 	# "Caudal_neurectoderm" = "#354E23",
# 	"Neurectoderm" = "#65A83E",
# 	"Neural_crest" = "#C3C388",
# 	"Forebrain_Midbrain_Hindbrain" = "#647a4f",
# 	"Spinal_cord" = "#CDE088",
# 	"Surface_ectoderm" = "#f7f79e",
# 	"Visceral_endoderm" = "#F6BFCB",
# 	"ExE_endoderm" = "#7F6874",
# 	"ExE_ectoderm" = "#989898",
# 	"Parietal_endoderm" = "#1A1A1A"
# )


opts$celltypes = c(
	"Epiblast",
	"Primitive_Streak",
	"Caudal_epiblast",
	"PGC",
	"Anterior_Primitive_Streak",
	"Notochord",
	"Def._endoderm",
	"Gut",
	"Nascent_mesoderm",
	"Mixed_mesoderm",
	"Intermediate_mesoderm",
	"Caudal_Mesoderm",
	"Paraxial_mesoderm",
	"Somitic_mesoderm",
	"Pharyngeal_mesoderm",
	"Cardiomyocytes",
	"Allantois",
	"ExE_mesoderm",
	"Mesenchyme",
	"Haematoendothelial_progenitors",
	"Endothelium",
	"Blood_progenitors_1",
	"Blood_progenitors_2",
	"Erythroid1",
	"Erythroid2",
	"Erythroid3",
	"NMP",
	"Rostral_neurectoderm",
	"Caudal_neurectoderm",
	"Neural_crest",
	"Forebrain_Midbrain_Hindbrain",
	"Spinal_cord",
	"Surface_ectoderm",
	"Visceral_endoderm",
	"ExE_endoderm",
	"ExE_ectoderm",
	"Parietal_endoderm"
)

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

opts$samples <- c(
  
  # E6.5
  "1", "5", "18",
  
  # E6.75
  "7",
  
  # "E7.0"
  "10", "14", "15", "30", "31", "32",
  
  # E7.25
  "23", "26", "27",
  
  # E7.5
  "2", "3", "4", "6", "19", "20",
  
  # E7.75
  "8", "9", "12", "13",
  
  # E8.0
  "16", "33", "34", "35",
  
  # E8.25
  "24", "25", "28",
  
  # E8.5
  "17", "29", "36", "37",
  
  # mixed_gastrulation
  "21", "22"
)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>%
  .[stripped==F & doublet==F] %>%
  .[,celltype:=factor(celltype, levels=names(opts$celltype.colors))]
  
