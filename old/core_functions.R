# scale_colour_annotation = function(...){
#   library(scales)
#   vals = c("PS/Mesendoderm" = "#efd5a0",#grey-brown ###
#            "ExEct 2" = "grey20",#darkgrey###
#            "Epiblast" = "#663300",#dark brown###
#            "Cardiac mesenchyme" = "thistle3",#light pink
#            "Neural crest" = "palegreen3",#light green
#            "Late parax/somit. meso." = "royalblue3",#blue
#            "Neuroectoderm" = "greenyellow",#midgreen
#            "Extraembryonic mesoderm" = "purple3",#purple
#            "Endothelium" = "orange",#orange
#            "Neural tube" = "olivedrab",#darkgreen
#            "ExEct 2" = "grey60",#light grey
#            "Hem-endo" = "firebrick1",#mid red
#            "NMP" = "#FAFF0A",#yellow
#            "Mixed mesoderm" = "navy",#navy
#            "Early parax. meso." = "steelblue3",#lightblue
#            "Cardiac" = "pink4",#dark pink
#            "ExEmb tissue" = "grey10",#???
#            "AVE/Def. endoderm" = "coral2",#dark goldenrod
#            "Prim. Endoderm" = "#A38566",#light brown ###
#            "ExEct doubs" = "grey40",#mid grey
#            "Mesoderm prog." = "#c4fffe",#skyblue
#            "Notochord"="turquoise",#bright blue
#            "Erythroid 1" = "firebrick3",#darkred
#            "Erythroid 2" = "red4")#scarlet
#   labels = names(vals)
#   names(vals) = NULL
#   
#   
#   
#   return(discrete_scale("fill", "Publication",
#                  manual_pal(vals),
#                  name = "All-cluster annotation",
#                  breaks = labels,
#                  labels = labels,
#                  drop = FALSE)
#   )
#   }
# 
# ggplot(mapping = aes(x = seq_along(all_colours), y = rep(1, length(all_colours)), fill = factor(seq_along(all_colours)))) +
#   geom_bar(stat = "identity") +
#   scale_colour_annotation() +
#   theme(legend.position = "right")

# CLUSTER TYPES
# CLUSTER TYPES
all_names = c("Epiblast",
              "PS/mesendoderm",
              "Erythroid 1",
              "NMPs",
              "Late neural tube/spinal cord",
              "ExE endoderm",
              "ExE mesoderm",
              "ExE ectoderm 1",
              "Neuroectoderm",
              "Hemato-endothelial",
              "Early parax. mesoderm",
              "ExE ectoderm 2",
              "Late parax. mesoderm",
              "AVE/def. endo/notochord",
              "Erythroid 2",
              "Late mixed mesoderm",
              "Visceral endoderm",
              "Early mixed mesoderm",
              "Parietal endoderm",
              "Neural crest/non-neural ectoderm")
names(all_names) = 1:length(all_names)

legend_order = match(c("Epiblast",
                       "PS/mesendoderm", 
                       "AVE/def. endo/notochord",
                       "Early mixed mesoderm",
                       "Early parax. mesoderm",
                       "Late parax. mesoderm",
                       "Late mixed mesoderm",
                       "Hemato-endothelial",
                       "Erythroid 1",
                       "Erythroid 2",
                       "ExE mesoderm",
                       "NMPs",
                       "Neuroectoderm",
                       "Neural crest/non-neural ectoderm",
                       "Late neural tube/spinal cord",
                       "ExE endoderm",
                       "Visceral endoderm",
                       "ExE ectoderm 1",
                       "ExE ectoderm 2",
                       "Parietal endoderm"
), all_names)

# COLOURS
all_colours = c("PS/mesendoderm" = "#efd5a0",#grey-brown ###
                "ExE ectoderm 2" = "grey20",#darkgrey###
                "Epiblast" = "#663300",#dark brown###
                "Neural crest/non-neural ectoderm" = "palegreen3",#light green
                "Late parax. mesoderm" = "royalblue3",#blue
                "Neuroectoderm" = "greenyellow",#midgreen
                "ExE mesoderm" = "purple3",#purple
                "Hemato-endothelial" = "orange",#orange
                "Late neural tube/spinal cord" = "olivedrab",#darkgreen
                "ExE ectoderm 1" = "grey60",#light grey
                "NMPs" = "#FAFF0A",#yellow
                "Late mixed mesoderm" = "navy",#navy
                "Early parax. mesoderm" = "steelblue1",#lightblue
                "Parietal endoderm" = "grey10",#very dark grey
                "AVE/def. endo/notochord" = "coral2",#dark goldenrod
                "ExE endoderm" = "plum4",#plum ###
                "Early mixed mesoderm" = "turquoise",#skyblue
                "Erythroid 1" = "firebrick3",#darkred
                "Erythroid 2" = "red4",
                "Visceral endoderm" = "lightpink1")#pink
all_colours = all_colours[match(all_names, names(all_colours))]
names(all_colours) = 1:length(all_colours)

celltype_colours = c(
  "Epiblast"	= "#683612",
  "PS/mesendoderm"	= "#DABE99",
  "Early mixed mesoderm"	= "#C594BF",
  "Early ExE mesoderm"	= "#DFCDE4",
  "ExE mesoderm"	= "#7253A2",
  "Allantois"	= "#532C8A",
  "Endothelium"	= "#B3793B",
  "Hemato-endothelial progenitors"	= "#FBBE92",
  "Erythroid 1"	= "#C72228",
  "Erythroid 2"	= "#EF4E22",
  "Cardiac mesenchyme"	= "#F7901D",
  "Cardiomyocytes"	= "#B51D8D",
  "Early paraxial mesoderm"	= "#3F84AA",
  "Late mixed mesoderm"	= "#C9EBFB",#pharyngeal mesoderm?
  "Intermediate mesoderm"	= "#139992",
  "Late parax. mesoderm"	= "#8DB5CE",
  "Somites"	= "#005579",
  "Early neurectoderm"	= "#A0CC47",
  "Forebrain"	= "#65A83E",
  "Midbrain/Hindbrain"	= "#354E23",
  "Pre-migratory neural crest"	= "#C3C388",#Cranial Neural Crest?
  "Neural crest"	= "#77783C",#Trunk Neural Crest?
  "Placodes"	= "#BBDCA8",
  "NMP"	= "#8EC792",
  "Spinal cord"	= "#CDE088",
  "Surface ectoderm"	= "#FFF574",
  "Notochord"	= "#0F4A9C",
  "PGC"	= "#FACB12",
  "Def. endoderm"	= "#F397C0",
  "Foregut"	= "#EF5A9D",
  "Midgut/Hindgut"	= "#CE4E82",
  "Visceral endoderm"	= "#F6BFCB",
  "Parietal endoderm"	= "#1A1A1A",
  "ExE endoderm"	= "#7F6874",
  "ExE ectoderm 1"	= "#989898",
  "ExE ectoderm 2"	= "#333333")


# ggplot(mapping = aes(x = seq_along(all_colours), y = rep(1, length(all_colours)), fill = factor(seq_along(all_colours), levels = legend_order))) +
#   geom_bar(stat = "identity") +
#   # scale_colour_annotation()
#   scale_fill_manual(values = all_colours, labels = all_names)
# 
# ggplot(data = plot_df, mapping = aes(x = V1, y = V2, col = factor(cluster))) +
#   geom_point(size = 1, alpha = 1) +
#   labs(x = "t-sne 1", y = "t-sne 2") +
#   scale_color_manual(values = all_colours, labels = all_names, drop = FALSE, name = "") +
#   theme_bw()


# ROUTINELY USED PACKAGES
library(irlba)

#load data to global variables
# sce_hvg: sce set object containing highly variable gene counts for all cells, excluding Xist and Y CHR genes
# genes: map to convert between Ensembl/MGI for the sce set above
# meta: latest metadata table
# sce_all: sce set containing counts for all genes
# genes_all: map to convert between Ensembl/MGI for the all-gene sceset
# mouse_ensembl: biomaRt Mart object for mmusculus
# gene_map: links genes (ensembl) to chromosomes.
load_data = function(normalise = TRUE, hvg = TRUE, corrected = FALSE, remove_doublets = FALSE){
  require(scran)
  require(scater)
  require(SingleCellExperiment)
  require(biomaRt)
  require(Matrix)
  
  counts = readRDS("/Users/ricard/data/gastrulation10x/data/raw_counts.rds")
  genes = read.table("/Users/ricard/data/gastrulation10x/data/genes.tsv", stringsAsFactors = F)
  meta <<- read.table("/Users/ricard/data/gastrulation10x/data/meta.tab", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  rownames(counts) = genes[,1] #ensembl
  colnames(counts) = meta$cell
  
  sce = SingleCellExperiment(assays = list("counts" = counts))
  
  if(normalise){
    sfs <<- read.table("/Users/ricard/data/gastrulation10x/data/sizefactors.tab", stringsAsFactors = F)[,1]
    sizeFactors(sce) = sfs
    sce = normalise(sce)
  }
  
  
  # sce_all <<- sce
  genes_all <<- genes
  
  if(remove_doublets){
    sce <<- normalise(sce_all[,!meta$doub.clust])
    meta <<- meta[!meta$doub.clust,]
  }
  
  if(hvg){
    # #remove Xist
    # sce = sce[-which(rownames(sce) == genes[match("Xist", genes[,2]),1])]
    # #remove Y
    # y_genes = gene_map$ensembl_gene_id[gene_map$chromosome_name == "Y"]
    # sce = sce[-which(rownames(sce) %in% y_genes),]
    
    # hvg.tab <<- read.table("/nfs/research1/marioni/jonny/embryos/data/hvgs.tab", stringsAsFactors = F, header = TRUE)
    # hvgs <<- hvg.tab[,1]
    # sce_hvg <<- sce[rownames(sce) %in% hvgs, ]
    # genes_hvg <<- genes_all[genes_all[,1] %in% hvgs, ]
    
    hvg.list <<- readRDS("/Users/ricard/data/gastrulation10x/data/hvg.list.rds")
  }
  
  if(corrected)
    sce_corrected <<- readRDS("/Users/ricard/data/gastrulation10x/data/sce_mnn.rds")
  
  return(sce)
}


#removed yellow from position 2 ("#FFFF00")
scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour", "Publication",
                 manual_pal(values = c(
                   "#000000", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                   "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                   "#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                   "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                   "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                   "#372101", "#FFB500", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                   "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                   "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                   "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                   "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                   "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                   "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                   "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")), ...)
}

#removed yellow from position 2 ("#FFFF00")
scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill", "Publication",
                 manual_pal(values = c(
                   "#000000", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                   "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                   "#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                   "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                   "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                   "#372101", "#FFB500", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                   "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                   "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                   "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                   "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                   "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                   "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                   "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")), ...)
}

smart_decomp = function(sce, stage, sample, loess.span = 0.05, assay = "logcounts"){
  require(scran)
  
  method = switch( (length(unique(sample)) > 1) + 1,
                  "base",
                  switch((length(unique(stage)) > 1) + 1,
                         "sample",
                         "stage_sample"))
  
  fit.base = trendVar(sce, use.spikes=FALSE, loess.args = list(span = loess.span), assay.type = assay)
  decomp.all = decomposeVar(sce, fit.base, assay.type = assay)

  if(method %in% c("sample", "stage_sample")){
    fit.sample <- trendVar(sce, use.spikes=FALSE, loess.args = list(span = loess.span), design = model.matrix(~ factor(sample)), assay.type = assay) 
    decomp.sample = decomposeVar(sce, fit.sample, assay.type = assay)
  }
  
  if(method == "stage_sample"){
    fit.stage <- trendVar(sce, use.spikes=FALSE, loess.args = list(span = loess.span), design = model.matrix(~ factor(stage)), assay.type = assay)
    decomp.stage = decomposeVar(sce, fit.stage, assay.type = assay)
  }
  
  
  if(method == "base"){
    decomp.new = decomp.all
  } else if(method == "sample"){
    decomp.new = decomp.sample
  } else {
    decomp.new = data.frame(mean = decomp.all$mean, 
                            total = decomp.all$total + decomp.sample$total - decomp.stage$total,
                            row.names = rownames(decomp.all),
                            stringsAsFactors = FALSE)
    
    #repeat the trend fitting like scran does
    #i.e. using log variance, so the loess can never fit a variance < 0
    fitvar = log(decomp.new$total[decomp.new$total>0])
    fitmean = decomp.new$mean[decomp.new$total>0]
    new_fit = loess(fitvar ~ fitmean, span = 0.05)
    decomp.new$tech = exp(predict(new_fit, data.frame(fitmean = decomp.new$mean)))
    decomp.new$tech[is.na(decomp.new$tech)] = 0
    
    #test via scran
    decomp.new$bio = decomp.new$total - decomp.new$tech
    mm = model.matrix(~stage + sample) #used both of these covariates
    decomp.new$p.value = testVar(decomp.new$total, decomp.new$tech, df = nrow(mm) - ncol(mm), test = "chisq")
    decomp.new$FDR = p.adjust(decomp.new$p.value, method = "fdr")
  }
  
  if(method %in% c("sample", "stage_sample")){
    #add the variance ratios to data frame
    #that is, the ratio of the variance of the mean between samples
    #to the mean of the variance within samples
    var_within = sapply(unique(sample), function(x){
      rowVars(as.matrix(logcounts(sce[,sample == x])))
    })
    means = sapply(unique(sample), function(x){
      Matrix::rowMeans(logcounts(sce[,sample == x]))
    })
    mean_var_within = rowMeans(var_within)
    var_mean = rowVars(means)
    decomp.new$ratio = rowVars(means)/rowMeans(var_within)
  }
  
  
  return(decomp.new)
}

choose_hvgs = function(decomp, fdr.pval = 0.01, nhvg = NULL, min.mean = 0.5e-3, max.mean = NULL, use.inflection = TRUE, return.plots=FALSE){
  require(biomaRt)
  require(ggplot2)
  require(gam)
  
  #retain for plot
  decomp.plot = as.data.frame(decomp)
  
  #remove unexpressed/lowly expressed
  decomp = decomp[decomp$mean > min.mean,]
  
  ratio_present = "ratio" %in% names(decomp)
  #override max.mean if inflection specified
  if(use.inflection & ratio_present){
    keep = !is.na(decomp$ratio)
    inflection_df = data.frame(gene = rownames(decomp)[keep], mean = decomp$mean[keep], ratio = decomp$ratio[keep])
    model_df = data.frame(logratio = log10(inflection_df$ratio),
                          logmean = log10(inflection_df$mean))
    #Fit a GAM model, see result in plots$gam_fit
    model = gam(data = model_df, formula = logratio ~ logmean + s(logmean, df = 10))
    #predict an evenly spread number of points across the data
    points = seq(from = min(log10(inflection_df$mean)),
                 to = max(log10(inflection_df$mean)),
                 length.out = 1e4)
    vals = predict(model, data.frame(logmean = points))
    #take the second derivative i.e. how is the gradient changing?
    twodiff = diff(diff(vals, 1), 1)
    names(twodiff) = points[-c(1:2)]
    #find the peak i.e. the steepest changing part on the peak at the end
    peak = which.max(twodiff)
    #low identifies when we start the ascent up to this peak
    #didn't use third derivative as it's unstable 
    # low = max(which(twodiff[seq_len(peak)]< 0)) #derivative hits 0
    low = max(which(twodiff[seq_len(peak)] < (max(twodiff) * 0.1) )) #derivative hits 10%
    #redefine max mean. [low] will be more conservative than [peak].
    max.mean = 10^(as.numeric(names(twodiff)[low]))
  } else if (use.inflection & !ratio_present){
    max.mean = as.numeric(quantile(decomp.plot$mean, 1-0.073)) #0.0932 from HVG script, mean value between stage hvg calls
  }
  
  #remove highly expressed genes
  if(!is.null(max.mean)){
    decomp = decomp[decomp$mean < max.mean,]
  }

  #remove sex
  db = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  gene_map = getBM(attributes = c("ensembl_gene_id", "chromosome_name"), filters = "ensembl_gene_id", values = rownames(decomp), mart = db)
  sex_genes = c("ENSMUSG00000086503", #Xist
                gene_map[gene_map[,2] == "Y",1])
  decomp = decomp[!rownames(decomp) %in% sex_genes,]
  
  #re-correct FDR (after thresholding)
  decomp$FDR = p.adjust(decomp$p.value, method = "fdr")
  
  #select HVGs
  if(all(is.null(fdr.pval), is.null(nhvg)) | (!is.null(fdr.pval) & !is.null(nhvg))){
    stop("Please use exactly one of fdr.pval or nhvg")
  }
  if(!is.null(fdr.pval)){
    decomp$select = decomp$FDR <= fdr.pval
  } else {
    decomp$select = FALSE
    decomp$select[order(decomp$p.value, -decomp$bio, decreasing = FALSE)[seq_len(nhvg)]] = TRUE
  }
  
  hvgs = rownames(decomp)[decomp$select]
  
  if(return.plots){
    decomp.plot$select = rownames(decomp.plot) %in% hvgs
    nhvg = sum(decomp.plot$select)
    
    chosen_genes = ggplot(decomp.plot, aes(x = mean, y = total)) +
      geom_point(aes(col = select), size= 0.8) +
      geom_line(aes(y = tech), col = "cornflowerblue") +
      scale_x_log10(breaks = c(1e-3, 1e-1, 1e1), labels = c("0.001", "0.1", "10")) + 
      scale_y_log10(breaks = c(1e-4, 1e-2, 1e0), labels = c("0.0001", "0.01", "1")) +
      theme_bw() +
      scale_color_manual(values = c("TRUE" = "coral", "FALSE" = "grey30")) +
      theme(legend.position = "none") +
      labs(x = "log10-count mean", y = "log10-count variance") +
      ggtitle(paste0("nHVG = ", nhvg))
    
    plots = list(chosen_genes = chosen_genes)
    
    if(use.inflection & ratio_present){
      
        gam_fit = ggplot(model_df, aes(x = logmean, y = logratio)) +
          geom_point(colour = "grey80") +
          geom_line(data = data.frame(X = points, Y = vals), mapping = aes(x = X, y = Y), col = "red") +
          theme_bw() +
          labs(x = "log10 mean(log-counts)", y = "log10 ratio") +
          geom_vline(mapping = aes(xintercept = log10(max.mean)), col = "coral", lty = "twodash")
        
        twodiff = ggplot(data.frame(x = points[-c(1:2)], y = twodiff), mapping = aes(x =x, y=y)) +
          geom_line() +
          labs(x = "log10 mean(log-counts)", y= "Second derivative") +
          theme_bw() +
          geom_vline(mapping = aes(xintercept = log10(max.mean)), col = "coral", lty = "twodash")
      
        plots = list(chosen_genes = chosen_genes, gam_fit = gam_fit, twodiff = twodiff)
    }
    return(list(hvgs = hvgs, plots = plots))
  }
  
  return(hvgs)
}


#from http://blog.schochastics.net/post/using-umap-in-r-with-rpython/
umap <- function(x,n_neighbors=10,min_dist=0.1,metric="euclidean"){
  x <- as.matrix(x)
  colnames(x) <- NULL
  rPython::python.exec( c( "def umap(data,n,mdist,metric):",
                           "\timport umap" ,
                           "\timport numpy",
                           "\tembedding = umap.UMAP(n_neighbors=n,min_dist=mdist,metric=metric).fit_transform(data)",
                           "\tres = embedding.tolist()",
                           "\treturn res"))
  
  res <- rPython::python.call( "umap", x,n_neighbors,min_dist,metric)
  do.call("rbind",res)
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}