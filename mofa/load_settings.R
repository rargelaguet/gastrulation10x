################
## Define I/O ##
################

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/gastrulation10x_mofa"
  io$gene.metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
} else if (grepl("ebi",Sys.info()['nodename'])) {
  io$basedir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x_mofa"
  io$gene.metadata <- "/hps/nobackup2/research/stegle/users/ricard/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
} else {
	stop("Computer not recognised")
}

io$pdfdir <- paste0(io$basedir,"/mofa/pdf")
io$sample.metadata <- paste0(io$basedir,"/sample_metadata.txt.gz")

####################
## Define options ##
####################

opts <- list()
opts$colors = c(
  "Epiblast" = "#635547",
	"Primitive Streak" = "#DABE99",
	"Caudal epiblast" = "#9e6762",
	
	"PGC" = "#FACB12",
	
	"Anterior Primitive Streak" = "#c19f70",
	"Notochord" = "#0F4A9C",
	"Def. endoderm" = "#F397C0",
	"Gut" = "#EF5A9D",
	
	"Nascent mesoderm" = "#C594BF",
	"Mixed mesoderm" = "#DFCDE4",
	"Intermediate mesoderm" = "#139992",
	"Caudal Mesoderm" = "#3F84AA",
	"Paraxial mesoderm" = "#8DB5CE",
	"Somitic mesoderm" = "#005579",
	"Pharyngeal mesoderm" = "#C9EBFB",
	"Cardiomyocytes" = "#B51D8D",
	"Allantois" = "#532C8A",
	"ExE mesoderm" = "#8870ad",
	"Mesenchyme" = "#cc7818",
	
	"Haematoendothelial progenitors" = "#FBBE92",
	"Endothelium" = "#ff891c",
	"Blood progenitors 1" = "#f9decf",
	"Blood progenitors 2" = "#c9a997",
	"Erythroid1" = "#C72228",
	"Erythroid2" = "#f79083",
	"Erythroid3" = "#EF4E22",
	
	"NMP" = "#8EC792",
	
	"Rostral neurectoderm" = "#65A83E",
	"Caudal neurectoderm" = "#354E23",
	"Neural crest" = "#C3C388",
	"Forebrain_Midbrain_Hindbrain" = "#647a4f",
	"Spinal cord" = "#CDE088",
	
	"Surface ectoderm" = "#f7f79e",
	
	"Visceral endoderm" = "#F6BFCB",
	"ExE endoderm" = "#7F6874",
	"ExE ectoderm" = "#989898",
	"Parietal endoderm" = "#1A1A1A"
)


opts$colors <- c(
  "Epiblast" = "grey70",
  "Primitive Streak" = "sandybrown",
  "Mesoderm" = "violetred",
  "ExE endoderm" = "#548B54",
  "ExE ectoderm" = "black"
)

# ignore...
# opts$factors.width <- 8
# opts$factors.height <- 4