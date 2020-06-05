wd <- "/Users/blancap/Documents/CAREER/GOTTGENS_partTime/Projects/BPS4_RArgelaguet/"

mat = read.csv(paste0(wd,"data/PAGA_connectivity.csv"),header=F)
labels = read.csv(paste0(wd,"data/PAGA_connectivity_labels.csv"))
rownames(mat)=colnames(mat)=labels$X0
matdist = 1-mat

diag(matdist) = 0
#print(lk)

write.csv(matdist,file=paste0(wd,"data/PAGA_distances.csv"))
#matdist = read.csv(paste0(wd,"data/PAGA_distances.csv"),header=T)

matdistX <- as.dist(matdist)

h <- hclust(matdistX, method = 'ward.D')
plot(h)


h2 <- hclust(matdistX, method = 'ward.D2')

plot(h2)


h2 <- hclust(matdistX, method = 'average')

plot(h2)
