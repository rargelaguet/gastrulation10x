library(data.table)
library(purrr)
library(ggplot2)

#####################
## Define settings ##
#####################

## Define I/O ##
io <- list()
io$stats <- "/Users/ricard/data/gastrulation10x/mofa/stats/r2.txt"
io$pdfdir <- "/Users/ricard/data/gastrulation10x/mofa/pdf"

## Define options ##
opts <- list()
opts$batch_size <- c( "0.25", "0.5" ) %>% as.numeric
opts$tau <- c( "0.05", "0.15", "0.25", "0.5", "0.75", "1.0" ) %>% as.numeric
opts$forgetting_rate <- c( "0.5", "0.25", "0" ) %>% as.numeric
opts$ntrials <- 1

#################
# Load results ##
#################

stats <- fread(io$stats) %>%
  # .[,tau:=as.factor(tau)] %>%
  # .[,forgetting_rate:=as.factor(forgetting_rate)] %>%
  # .[,time:=time/60] %>% .[time<600] %>%
  # .[,batch_size:=as.numeric(batch_size)] %>%
  .[trial==1] %>% .[,trial:=as.factor(trial)]

stats_nostochastic <- stats[batch_size=="nostochastic"]
stats_stochastic <- stats[tau%in%opts$tau & forgetting_rate%in%opts$forgetting_rate & batch_size%in%opts$batch_size] %>%
  .[,batch_size:=as.factor(sprintf("Batch size: %d%%",round(100*as.numeric(batch_size))))] %>%
  .[,tau:=sprintf("Tau: %s",tau)] 



##########
## Plot ##
##########

tmp <- stats_stochastic %>%
  .[,batchsize_tau:=paste(tau,batch_size,sep="   ")]

p <- ggplot(tmp, aes(x=as.factor(forgetting_rate), y=r2)) +
  geom_bar(stat="identity", position="dodge", fill="darkgrey", color="black", alpha=0.75) +
  geom_hline(yintercept=stats_nostochastic$r2, linetype="dashed") +
  # geom_hline(yintercept=0.80, linetype="dashed", size=0.8) +
  facet_wrap(~batchsize_tau, scales = "free_x", ncol=2) +
  scale_y_continuous(limits=c(0,1)) +
  labs(x="Forgetting rate", y="Variance explained (R2)") +
  theme(
    panel.spacing = unit(1.5, "lines"),
    axis.title.x = element_text(colour="black", size=15),
    axis.title.y = element_text(colour="black", size=15),
    axis.text.x = element_text(colour="black",size=rel(1.1)),
    axis.text.y = element_text(colour="black",size=rel(1.1)),
    axis.ticks = element_line(colour="black"),
    axis.line = element_line(color="black"),
    legend.position="none",
    strip.text = element_text(size=rel(1.0), color="black"),
    # legend.title = element_blank(),
    legend.direction = "vertical",
    legend.key.width=unit(1.2,"line"),
    legend.key.height=unit(1.0,"line"),
    legend.text = element_text(size=15),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )
print(p)

pdf(paste0(io$pdfdir,"/r2.pdf"), useDingbats = F, onefile = F, width=7.5, height=9)
print(p)
dev.off()