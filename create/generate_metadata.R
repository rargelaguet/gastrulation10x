meta <- fread("/Users/ricard/data/gastrulation10x/meta.tab") %>%
  .[,c("cell","barcode","sample","stage","theiler","doub.density","doublet","stripped","celltype","colour")]

fwrite(meta, "/Users/ricard/data/gastrulation10x/sample_metadata.txt", quote=F, sep="\t", col.names=T)
colnames(meta)
