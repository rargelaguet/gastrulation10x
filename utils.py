import anndata as anndata
import scanpy as sc

def load_adata(adata_file, metadata_file = None, normalise = False, cells = None, features = None, filter_lowly_expressed_genes = False, set_colors = True):

	adata = sc.read(adata_file)

	if metadata_file is not None:
		metadata = pd.read_table(metadata_file, delimiter="\t", header=0)
		assert np.all(metadata.cell.isin(adata.obs.index))
		assert metadata.shape[0] == adata.shape[0]
		adata.obs = metadata.set_index("cell").reindex(adata.obs.index)

	# metadata_to_anndata = (adata.obs >> 
	#     select("barcode","sample") >>
	#     left_join(metadata, by=["barcode","sample"])
	# )
	# metadata_to_anndata = metadata_to_anndata.set_index("cell")
	# assert metadata_to_anndata.shape[0] == adata.shape[0]
	# metadata_to_anndata.head()
	# adata.obs = metadata_to_anndata

	if cells is not None:
		adata = adata[cells,:]

	if features is not None:
		adata = adata[:,features]

	if filter_lowly_expressed_genes:
		sc.pp.filter_genes(adata, min_counts=10)

	if normalise:
		sc.pp.normalize_total(adata, target_sum=None, exclude_highly_expressed=False)
		sc.pp.log1p(adata)

	if set_colors:
		colPalette_celltypes = [opts["celltype_colors"][i.replace(" ","_").replace("/","_")] for i in sorted(np.unique(adata.obs['celltype']))]
		adata.uns['celltype_colors'] = colPalette_celltypes
		colPalette_stages = [opts["stages_colors"][i.replace(" ","_").replace("/","_")] for i in sorted(np.unique(adata.obs['stage']))]
		adata.uns['stage_colors'] = colPalette_stages

	return adata

