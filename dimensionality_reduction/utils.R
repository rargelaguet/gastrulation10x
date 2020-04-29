# aggregated_celltype_colours <- c(
#   "Epiblast"="gray70",
#   "Mesoderm"="#CD3278",
#   "Primitive Streak"="sandybrown",
#   "Endoderm"="#43CD80",
#   "Ectoderm"="steelblue",
#   "ExE Ectoderm"="black",
#   "Blood"="darkred"
# )

# celltype_colours = c(
#  "Epiblast" = "#635547",
#  "Primitive Streak" = "#DABE99",
#  "Caudal epiblast" = "#9e6762",
#  "PGC" = "#FACB12",
#  "Anterior Primitive Streak" = "#c19f70",
#  "Notochord" = "#0F4A9C",
#  "Def. endoderm" = "#F397C0",
#  "Gut" = "#EF5A9D",
#  "Nascent mesoderm" = "#C594BF",
#  "Mixed mesoderm" = "#DFCDE4",
#  "Intermediate mesoderm" = "#139992",
#  "Caudal Mesoderm" = "#3F84AA",
#  "Paraxial mesoderm" = "#8DB5CE",
#  "Somitic mesoderm" = "#005579",
#  "Pharyngeal mesoderm" = "#C9EBFB",
#  "Cardiomyocytes" = "#B51D8D",
#  "Allantois" = "#532C8A",
#  "ExE mesoderm" = "#8870ad",
#  "Mesenchyme" = "#cc7818",
#  "Haematoendothelial progenitors" = "#FBBE92",
#  "Endothelium" = "#ff891c",
#  "Blood progenitors 1" = "#f9decf",
#  "Blood progenitors 2" = "#c9a997",
#  "Erythroid1" = "#C72228",
#  "Erythroid2" = "#f79083",
#  "Erythroid3" = "#EF4E22",
#  "NMP" = "#8EC792",
#  "Rostral neurectoderm" = "#65A83E",
#  "Caudal neurectoderm" = "#354E23",
#  "Neural crest" = "#C3C388",
#  "Forebrain_Midbrain_Hindbrain" = "#647a4f",
#  "Spinal cord" = "#CDE088",
#  "Surface ectoderm" = "#f7f79e",
#  "Visceral endoderm" = "#F6BFCB",
#  "ExE endoderm" = "#7F6874",
#  "ExE ectoderm" = "#989898",
#  "Parietal endoderm" = "#1A1A1A"
# )

# aggregated_celltypes = c(
#  "Epiblast" = "Epiblast",
#  "Primitive Streak" = "Primitive Streak",
#  "Caudal epiblast" = "Primitive Streak",
#  "PGC" = "PGC",
#  "Anterior Primitive Streak" = "Primitive Streak",
#  "Notochord" = "Endoderm",
#  "Def. endoderm" = "Endoderm",
#  "Gut" = "Endoderm",
#  "Nascent mesoderm" = "Mesoderm",
#  "Mixed mesoderm" = "Mesoderm",
#  "Intermediate mesoderm" = "Mesoderm",
#  "Caudal Mesoderm" = "Mesoderm",
#  "Paraxial mesoderm" = "Mesoderm",
#  "Somitic mesoderm" = "Mesoderm",
#  "Pharyngeal mesoderm" = "Mesoderm",
#  "Cardiomyocytes" = "Mesoderm",
#  "Allantois" = "Mesoderm",
#  "ExE mesoderm" = "Mesoderm",
#  "Mesenchyme" = "Mesoderm",
#  "Haematoendothelial progenitors" = "Blood",
#  "Endothelium" = "Blood",
#  "Blood progenitors 1" = "Blood",
#  "Blood progenitors 2" = "Blood",
#  "Erythroid1" = "Blood",
#  "Erythroid2" = "Blood",
#  "Erythroid3" = "Blood",
#  "NMP" = "Ectoderm",
#  "Rostral neurectoderm" = "Ectoderm",
#  "Caudal neurectoderm" = "Ectoderm",
#  "Neural crest" = "Ectoderm",
#  "Forebrain_Midbrain_Hindbrain" = "Ectoderm",
#  "Spinal cord" = "Ectoderm",
#  "Surface ectoderm" = "Ectoderm",
#  "Visceral endoderm" = "Endoderm",
#  "ExE endoderm" = "Endoderm",
#  "ExE ectoderm" = "ExE Ectoderm",
#  "Parietal endoderm" = "Endoderm"
# )

aggregated_celltypes = c(
 "Epiblast" = "Epiblast",
 "Primitive Streak" = "Primitive Streak",
 "Caudal epiblast" = "Primitive Streak",
 "PGC" = "PGC",
 "Anterior Primitive Streak" = "Primitive Streak",
 "Notochord" = "Notochord",
 "Def. endoderm" = "Endoderm",
 "Gut" = "Endoderm",
 "Nascent mesoderm" = "Nascent mesoderm",
 "Mixed mesoderm" = "Mesoderm",
 "Intermediate mesoderm" = "Mesoderm",
 "Caudal Mesoderm" = "Mesoderm",
 "Paraxial mesoderm" = "Mesoderm",
 "Somitic mesoderm" = "Mesoderm",
 "Pharyngeal mesoderm" = "Mesoderm",
 "Cardiomyocytes" = "Cardiomyocytes",
 "Allantois" = "Allantois",
 "ExE mesoderm" = "ExE mesoderm",
 "Mesenchyme" = "Mesenchyme",
 "Haematoendothelial progenitors" = "Haematoendothelial progenitors",
 "Endothelium" = "Endothelium",
 "Blood progenitors 1" = "Blood progenitors",
 "Blood progenitors 2" = "Blood progenitors",
 "Erythroid1" = "Erythroid",
 "Erythroid2" = "Erythroid",
 "Erythroid3" = "Erythroid",
 "NMP" = "NMP",
 "Rostral neurectoderm" = "Ectoderm",
 "Caudal neurectoderm" = "Ectoderm",
 "Neural crest" = "Neural crest",
 "Forebrain_Midbrain_Hindbrain" = "Forebrain_Midbrain_Hindbrain",
 "Spinal cord" = "Spinal cord",
 "Surface ectoderm" = "Ectoderm",
 "Visceral endoderm" = "ExE endoderm",
 "ExE endoderm" = "ExE endoderm",
 "ExE ectoderm" = "ExE Ectoderm",
 "Parietal endoderm" = "ExE endoderm"
)

colors <- c(
  "Ectoderm" = "steelblue",
  "Mesoderm" = "#CD3278",
  "Endoderm" = "#43CD80",
  "Blood" = "#C72228",

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
  "Blood progenitors" = "#f9decf",
  "Blood progenitors 1" = "#f9decf",
  "Blood progenitors 2" = "#c9a997",
  "Erythroid1" = "#C72228",
  "Erythroid2" = "#f79083",
  "Erythroid" = "#f79083",
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
plot.dimred <- function(plot_df) {
  ggplot(data=plot_df, mapping=aes(x=V1, y=V2)) +
    ggrastr::geom_point_rast(aes(size=mapped, alpha=mapped, colour=mapped)) +
    scale_size_manual(values = c("TRUE"=opts$size.mapped, "FALSE"=opts$size.nomapped)) +
    scale_alpha_manual(values = c("TRUE"=opts$alpha.mapped, "FALSE"=opts$alpha.nomapped)) +
    scale_colour_manual(values = c("TRUE"="red", "FALSE"="grey")) +
    guides(colour = guide_legend(override.aes = list(size=6))) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )
}