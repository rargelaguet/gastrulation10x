############
## Ribbon ##
############

#Each of these is a vector, with each entry corresponding to a single cell.
#I.e. row one of clusters and timepoints corresponds to the SAME cell.
#I.e., use metadata columns!
#timepoints are x-axis, ordered as you want them to appear
#clusters are y-axis shading 
#colours should be a named vector of colours, names are the clusters specified.
#the order of colours specifies the order in the plot
plot_fn = function(timepoints, clusters, colours = celltype_colours, plateau = FALSE, embryo_doubling = FALSE){
  
  df = data.frame(
    cluster = rep(levels(clusters), length(levels(timepoints))),
    stage = do.call(c, lapply(as.character(levels(timepoints)), rep, times = length(levels(clusters))))
  )
  
  df$ranking = match(df$cluster, names(colours))
  df= df[order(df$stage, df$ranking),]
  df$frac = sapply(seq_len(nrow(df)), function(x){
    return(sum(clusters == df$cluster[x] & timepoints == df$stage[x])/sum(timepoints == df$stage[x]))
  })
  
  df$cumfrac = NA
  for(x in 1:nrow(df)){
    df$cumfrac[x] = sum(df$frac[df$stage == df$stage[x] & df$ranking < df$ranking[x]])
  }
  
  df$xpos = match(df$stage, levels(timepoints))
  if(plateau){
    df1 = df
    df2 = df
    df1$xpos = df1$xpos - 0.2
    df2$xpos = df2$xpos + 0.2
    df = rbind(df1, df2)
  }
  
  p = ggplot(df, aes(x = xpos, 
                     # x = factor(stage, levels = levels(stage)[order(levels(stage))]), 
                     ymin = cumfrac, 
                     ymax = cumfrac + frac, 
                     fill = factor(cluster, levels = names(colours)), 
                     col = factor(cluster, levels = names(colours)))) +
    geom_ribbon() +
    scale_fill_manual(values = colours, drop = FALSE, name = "") +
    scale_color_manual(values = colours, drop = FALSE, name = "") +
    scale_x_continuous(breaks = seq_along(levels(timepoints)), labels = levels(timepoints), name = "")
  
  return(p)
  
}




source("/Users/ricard/gastrulation10x/settings.R")

sample_metadata <- sample_metadata[!stage=="mixed_gastrulation"]

timepoints <- factor(sample_metadata$stage, levels=opts$stages[opts$stages!="mixed_gastrulation"])
clusters <- sample_metadata$celltype
colours <- opts$celltype.colors

plot_fn(timepoints, clusters, colours, plateau = TRUE) +
  theme_classic() +
  theme(
    # axis.text.x = element_text(color="black"),
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )
