##======  Read graph abstraction coordinates ======##


coords_file = wd + "scripts/figures/graph_abstraction/paga_graph_coords_v6_20180913_mod.gdf"

coords = pd.read_csv(coords_file, sep=",",names=["name_VARCHAR",
                                               
                                                 "label_VARCHAR",
                                                                      'height_DOUBLE',
                                                                      'x_DOUBLE',
                                                                     'y_DOUBLE',
                                                                      'color_VARCHAR1',
                                                                      'col2',
                                                                      'col3','col4'])

coords_sorted = coords


coords_x = [float(i) for i in coords_sorted.loc[:,'x_DOUBLE']]
#print(graph_x)
coords_y = [float(i) for i in coords_sorted.loc[:,'y_DOUBLE']]
#print(graph_y)

coords_matrix = np.array([coords_x, coords_y])
print(coords_matrix.T.shape)
#Add coordinates in to adata file



##======  Plot PAGA coloured by celltypes subclusters ======##

sc.pl.paga(adata, threshold=0.23,fontsize=3,edge_width_scale=0.15,node_size_power=0.5,color='celltype_new',
           #save="clustersMergedCelltypes_v2_AGA.svg",export_to_gexf=True,
           pos=coords_matrix.T)


connect = adata.uns['paga']['connectivities'].todense()
np.savetxt(wd + 'scripts/write/BPS4_RArguelaguet/PAGA_connectivity.csv', connect, delimiter=',')

d = pd.DataFrame(data=np.unique(adata.obs['celltype_new']))
d.to_csv(wd + 'scripts/write/BPS4_RArguelaguet/PAGA_connectivity_labels.csv', sep=',')


##======  Plot PAGA coloured by stages ======##

colStage = [spectralDic[str(i)] for i in list(adata.obs['stage'])]
adata.uns['stage_colors'] = spectralPal
stageFollow = np.unique(adata.obs['stage'])

spectralPalFollow = [spectralDic[i] for i in stageFollow]


adata.uns['stage_colors'] = spectralPalFollow
paga2(adata,
           threshold=0.23,color='stage',fontsize=1,edge_width_scale=0.15,node_size_power=0.25,
           #save="clustersMergedCelltypes_v2_stage_AGA.svg",
      pos=coords_matrix.T
          )
sc.pl.paga(adata, threshold=0.23,fontsize=3,edge_width_scale=0.15,node_size_power=0.5,
           color=['ENSMUSG00000012396'],cmap=cmap,
           save="clustersMergedCelltypes_v2_AGA_Nanog.svg",
           export_to_gexf=True,
           pos=coords_matrix.T)

sc.pl.paga(adata, threshold=0.23,fontsize=3,edge_width_scale=0.15,node_size_power=0.5,
           color=['ENSMUSG00000012396'],cmap=cmap,
           save="clustersMergedCelltypes_v2_AGA_Nanog.pdf",
           export_to_gexf=True,
           pos=coords_matrix.T)


sc.pl.paga(adata, threshold=0.23,fontsize=3,edge_width_scale=0.15,node_size_power=0.5,
           color=['ENSMUSG00000037025'],cmap=cmap,
           save="clustersMergedCelltypes_v2_AGA_Foxa2.svg",
           export_to_gexf=True,
           pos=coords_matrix.T)

sc.pl.paga(adata, threshold=0.23,fontsize=3,edge_width_scale=0.15,node_size_power=0.5,
           color=['ENSMUSG00000037025'],cmap=cmap,
           save="clustersMergedCelltypes_v2_AGA_Foxa2.pdf",
           export_to_gexf=True,
           pos=coords_matrix.T)

sc.pl.paga(adata, threshold=0.23,fontsize=3,edge_width_scale=0.15,node_size_power=0.5,
           color=['ENSMUSG00000025219'],cmap=cmap,
           save="clustersMergedCelltypes_v2_AGA_Fgf8.svg",
           export_to_gexf=True,
           pos=coords_matrix.T)

sc.pl.paga(adata, threshold=0.23,fontsize=3,edge_width_scale=0.15,node_size_power=0.5,
           color=['ENSMUSG00000025219'],cmap=cmap,
           save="clustersMergedCelltypes_v2_AGA_Fgf8.pdf",
           export_to_gexf=True,
           pos=coords_matrix.T)