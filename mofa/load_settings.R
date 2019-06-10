## Define I/O ##
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/gastrulation10x"
  io$gene.metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
} else {
  io$basedir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x"
  io$gene.metadata <- "/hps/nobackup2/research/stegle/users/ricard/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
}
io$pdfdir <- paste0(io$basedir,"/mofa/pdf")
io$sample.metadata <- paste0(io$basedir,"/sample_metadata.txt")

## Define options ##
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
	"Forebrain/Midbrain/Hindbrain" = "#647a4f",
	"Spinal cord" = "#CDE088",
	
	"Surface ectoderm" = "#f7f79e",
	
	"Visceral endoderm" = "#F6BFCB",
	"ExE endoderm" = "#7F6874",
	"ExE ectoderm" = "#989898",
	"Parietal endoderm" = "#1A1A1A"
)

# counts = logcounts
# sample_order should just be the sample numbers from biggest to smallest
# and samples and timepoints the metadata columns
# doBatchCorrect = function(
#   counts, timepoints, samples, 
#   timepoint_order = c("E6.5", "E6.75", "E7.0", "mixed_gastrulation", "E7.25", "E7.5", "E7.75", "E8.0", "E8.25", "E8.5"), 
#   sample_order, 
#   npc = 50,
#   pc_override = NULL, 
#   BPPARAM = SerialParam()) {
#   require(scran)
#   require(irlba)
#   require(BiocParallel)
#   
#   if(!is.null(pc_override)){
#     pca = pc_override
#   } else {
#     pca = prcomp_irlba(t(counts), n = npc)$x
#     rownames(pca) = colnames(counts)
#   }
#   
#   if(length(unique(samples)) == 1){
#     return(pca)
#   }
#   
#   #create nested list
#   pc_list = lapply(unique(timepoints), function(tp){
#     sub_pc = pca[timepoints == tp, , drop = FALSE]
#     sub_samp = samples[timepoints == tp]
#     list = lapply(unique(sub_samp), function(samp){
#       sub_pc[sub_samp == samp, , drop = FALSE]
#     })
#     names(list) = unique(sub_samp)
#     return(list)
#   })
#   
#   names(pc_list) = unique(timepoints)
#   
#   # arrange to match timepoint order
#   pc_list = pc_list[order(match(names(pc_list), timepoint_order))]
#   pc_list = lapply(pc_list, function(x){
#     x[order(match(names(x), sample_order))]
#   })
#   
#   # perform corrections within list elements (i.e. within stages)
#   correct_list = lapply(pc_list, function(x){
#     if(length(x) > 1){
#       return(do.call(fastMNN, c(x, "pc.input" = TRUE, BPPARAM = BPPARAM))$corrected)
#     } else {
#       return(x[[1]])
#     }
#   })
#   
#   # perform correction over list
#   if(length(correct_list)>1){
#     correct = do.call(fastMNN, c(correct_list, "pc.input" = TRUE, BPPARAM = BPPARAM))$corrected
#   } else {
#     correct = correct_list[[1]]
#   }
#   
#   correct = correct[match(colnames(counts), rownames(correct)),]
#   
#   return(correct)
#   
# }