# import numpy as np
import anndata as anndata
# import scipy as s
import scanpy as sc
# import matplotlib.pyplot as pl
# import seaborn
import pandas as pd
# import csv
import os
import dfply
from re import search

#########
## I/O ##
#########

io = {}

host = os.uname()[1]
if search("ricard", host):
	io["basedir"] = "/Users/ricard/data/gastrulation10x"
	io["gene_metadata"] = "/Users/ricard/data/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt"
elif search("ebi", host):
	io["basedir"] = "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x"
	io["gene_metadata"] = "/hps/nobackup2/research/stegle/users/ricard/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
else:
	print("Computer not recognised"); exit()

io["metadata"] = io["basedir"] + "/sample_metadata.txt.gz"

# doublets and poor quality cells already removed
io["anndata"] = io["basedir"] + "/processed/scanpy/AnnData.h5ad"

#############
## Options ##
#############

opts = {}

opts["samples"] = [

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
]

opts["stages"] = [
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
]


opts["celltypes"] = [
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
]

opts["celltype_colors"] = {
	"Epiblast" : "#635547",
	"Primitive_Streak" : "#DABE99",
	"Caudal_epiblast" : "#9e6762",
	"PGC" : "#FACB12",
	"Anterior_Primitive_Streak" : "#c19f70",
	"Notochord" : "#0F4A9C",
	"Def._endoderm" : "#F397C0",
	"Gut" : "#EF5A9D",
	"Nascent_mesoderm" : "#C594BF",
	"Mixed_mesoderm" : "#DFCDE4",
	"Intermediate_mesoderm" : "#139992",
	"Caudal_Mesoderm" : "#3F84AA",
	"Paraxial_mesoderm" : "#8DB5CE",
	"Somitic_mesoderm" : "#005579",
	"Pharyngeal_mesoderm" : "#C9EBFB",
	"Cardiomyocytes" : "#B51D8D",
	"Allantois" : "#532C8A",
	"ExE_mesoderm" : "#8870ad",
	"Mesenchyme" : "#cc7818",
	"Haematoendothelial_progenitors" : "#FBBE92",
	"Endothelium" : "#ff891c",
	"Blood_progenitors" : "#c9a997",
	"Blood_progenitors_1" : "#f9decf",
	"Blood_progenitors_2" : "#c9a997",
	"Erythroid" : "#EF4E22",
	"Erythroid1" : "#C72228",
	"Erythroid2" : "#f79083",
	"Erythroid3" : "#EF4E22",
	"NMP" : "#8EC792",
	"Neurectoderm" : "#65A83E",
	"Rostral_neurectoderm" : "#65A83E",
	"Caudal_neurectoderm" : "#354E23",
	"Neural_crest" : "#C3C388",
	"Forebrain_Midbrain_Hindbrain" : "#647a4f",
	"Spinal_cord" : "#CDE088",
	"Surface_ectoderm" : "#f7f79e",
	"Visceral_endoderm" : "#F6BFCB",
	"ExE_endoderm" : "#7F6874",
	"ExE_ectoderm" : "#989898",
	"Parietal_endoderm" : "#1A1A1A"
}

opts["stages_colors"] = {
	'E6.5':"#D53E4F",
	'E6.75':"#F46D43",
	'E7.0':"#FDAE61",
	'E7.5':"#FFFFBF",
	'E7.25':"#FEE08B",
	'E7.75':"#E6F598",
	'E8.0':"#ABDDA4",
	'E8.5':"#3288BD",
	'E8.25':"#66C2A5",
	'mixed_gastrulation': "#A9A9A9"  
	
}

###############
## Load data ##
###############

# Load AnnData
# adata = sc.read(io["anndata"])

# Load sample metadata
# metadata = pd.read_table(io["metadata"], delimiter="\t", header=0)

# Add metadata to AnnData object
# assert metadata.shape[0] == adata.shape[0]
# metadata = metadata.set_index("cell").loc[adata.obs.index,:]
# adata.obs = metadata

################
## Set colors ##
################

# colPalette_celltypes = [opts["celltype_colors"][i.replace(" ","_")] for i in sorted(np.unique(adata.obs['celltype']))]
# adata.uns['celltype_colors'] = colPalette_celltypes

# colPalette_stages = [opts["stages_colors"][i.replace(" ","_")] for i in sorted(np.unique(adata.obs['stage']))]
# adata.uns['stage_colors'] = colPalette_stages

