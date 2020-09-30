
#install.packages("bnlearn")
#install.packages("BiocManager")
#BiocManager::install("Rgraphviz")

library(bnlearn)
library(parallel)
library(Rgraphviz)

#
## Wallace 
#
counts <- read.csv(file="dicho_Shi_samples_genusmicro_aclame_iceberg_megaresgroup.csv", row.names=1,header=TRUE)

columns <- colnames(counts)
counts[, columns] <- lapply(columns, function(x) as.numeric(counts[[x]]))
counts.dedup <- dedup(counts, .95, debug = FALSE)

set.seed(42)
sink(file = "learning-log-hc.boot.strength.txt")

dag.hybrid.group <- hc(counts)
cl <- parallel::makePSOCKcluster(16,outfile="debug.txt")
boot.hc.hybrid.group = boot.strength(data=counts,R=100,algorithm="hc", algorithm.args = list(cluster=cl))

sink()

stopCluster(cl)

plot(boot.hc.hybrid.group)

avg.boot.hc.hybrid.group <- averaged.network(boot.hc.hybrid.group, threshold = .7)
avg.boot.hc.dag.hybrid.group <- cextend(avg.boot.hc.hybrid.group)
#pdag2dag(avg.boot.hc.dag.hybrid.group)
fitted.hybrid.group <- bn.fit(avg.boot.hc.dag.hybrid.group,counts)

#Visualize it:
# d.gph.hybrid.group <- graphviz.plot(avg.boot.hc.dag.hybrid.group)
# pdf(file="final_network_hybrid_group.pdf",width=20,height=20)
# plot(d.gph.hybrid.group, nodeAttrs=makeNodeAttrs(d.gph.hybrid.group, fontsize=34))
# dev.off()

save.image(file="network_Shi_100rep_dicho_all_rumen_data.RData")
