library(data.table)
library(purrr)
library(ggplot2)

## Define I/O ##
io <- list()
io$stats <- "/Users/ricard/data/gastrulation10x/mofa/stats/stats.txt.gz"
io$pdfdir <- "/Users/ricard/data/gastrulation10x/mofa/pdf"

## Define options ##
opts <- list()
# opts$batch_size <- c( 0.25, 0.5, 1.0 )
# opts$tau <- c( 0.05, 0.15, 0.25, 0.5, 0.75, 1.0 )
# opts$forgetting_rate <- c( 0.5, 0.25, 0 )
opts$batch_size <- c( "0.25", "0.5" )
opts$tau <- c( "0.05", "0.15", "0.25", "0.5", "0.75", "1.0" )
opts$forgetting_rate <- c( "0.5", "0.25", "0" )
opts$ntrials <- 1


stats <- fread(sprintf(cmd="zcat < %s",io$stats)) %>%
  # .[,tau:=as.factor(tau)] %>%
  # .[,forgetting_rate:=as.factor(forgetting_rate)] %>%
  # .[,batch_size:=as.numeric(batch_size)] %>%
  .[trial==1] %>% #.[,trial:=as.factor(trial)] %>%
  .[complete.cases(.)]

stats_nostochastic <- stats[batch_size=="nostochastic"]
stats_stochastic <- stats[tau%in%opts$tau & forgetting_rate%in%opts$forgetting_rate & batch_size%in%opts$batch_size] %>%
  .[,batch_size:=as.factor(sprintf("Batch size: %d%%",round(100*as.numeric(batch_size))))] %>%
  .[,tau:=sprintf("Tau: %s",tau)] 


tmp <- stats_stochastic# %>% .[iter>200]

tmp[,batchsize_tau:=paste(tau,batch_size,sep="   ")]

# medium_convergence <- 0.00001
# slow_convergence <- 0.000001

tmp[,deltaELBO:=abs(100*((elbo-data.table::shift(elbo, 1, type="lag", fill=NA))/elbo[1])), by=c("batchsize_tau","trial")]
# foo <- tmp[deltaELBO<fast_convergence, .SD[which.min(iter)], c("batchsize_tau","trial")]
# bar <- tmp[deltaELBO<medium_convergence, .SD[which.min(iter)], c("batchsize_tau","trial")]
# baz <- tmp[deltaELBO<slow_convergence, .SD[which.min(iter)], c("batchsize_tau","trial")]

tmp <- tmp[iter!=1]

p <- ggplot(tmp, aes(x=iter, y=log(elbo), color=as.factor(forgetting_rate))) +
  # geom_bar(stat="identity", position="dodge") +
  geom_line() +
  geom_hline(yintercept=log(max(stats_nostochastic$elbo)), linetype="dashed") +
  # geom_vline(aes(xintercept=iter), linetype="dashed", data=baz, color="grey") +
  facet_wrap(~batchsize_tau, scales = "free_x", ncol=2) +
  # scale_y_continuous(expand=c(0,0)) +
  labs(x="Iteration", y="log Evidence Lower Bound") +
  theme(
    panel.spacing = unit(1.5, "lines"),
    axis.title.x = element_text(colour="black", size=15),
    axis.title.y = element_text(colour="black", size=15),
    axis.text.x = element_text(colour="black",size=rel(1.1)),
    axis.text.y = element_text(colour="black",size=rel(1.1)),
    axis.ticks = element_line(colour="black"),
    axis.line = element_line(color="black"),
    legend.position="right",
    strip.text = element_text(size=rel(1.0), color="black"),
    # legend.title = element_blank(),
    legend.direction = "vertical",
    legend.key.width=unit(1.2,"line"),
    legend.key.height=unit(1.0,"line"),
    legend.text = element_text(size=rel(1.1)),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )
print(p)

pdf(paste0(io$pdfdir,"/elbo.pdf"), useDingbats = F, onefile = F, width=16, height=17)
print(p)
dev.off()