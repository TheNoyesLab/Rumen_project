## Script to summarize counts
### Code for dataset analysis 
library(data.table)
library(tidyr)
#setwd("~/Dropbox/WRITING/NCBA2_JAN2018/NCBA2_analysis/") 
#source('scripts/Frontiers_NCBA2_analysis.R')
options(scipen = 999) # to decrease the use of scientific notation



setkey(amr_melted_raw_analytic,ID) 
setkey(amr_melted_analytic,ID) 
setkey(ACLAME_melted_analytic,ID)

setkey(metadata,ID)

amr_melted_analytic <- amr_melted_analytic[metadata]
amr_melted_raw_analytic <- amr_melted_raw_analytic[metadata]

ACLAME_melted_analytic2 <- ACLAME_melted_analytic[metadata]

ACLAME_melted_analytic


##### AMR exploratory #####

## By ID, faceted by Study, Class level
AMR_class_sum <- ACLAME_melted_analytic2[Level_ID=="Gene", .(sum_class= sum(Normalized_Count), Species),by=.(ID, Name, Study)]
AMR_class_sum[,total:= sum(sum_class), by=.(ID)]
AMR_class_sum[,percentage:= sum_class/total ,by=.(ID, Name) ]
# only keep taxa greater than 1%
AMR_class_sum <- AMR_class_sum[percentage > .1]
AMR_class_sum[,total:= sum(sum_class), by=.(ID)]
AMR_class_sum[,percentage:= sum_class/total ,by=.(ID, Name) ]




AMR_class_sum$Gene <- AMR_class_sum$Name
#AMR_class_sum[,percentage:= round(sum_class/total, digits=2) ,by=.(ID, Name) ] removes some with low proportions
ggplot(AMR_class_sum, aes(x = ID, y = percentage, fill = Gene)) + 
  geom_bar(stat = "identity")+
  facet_wrap( ~ Study, scales='free',ncol = 2) +
  #scale_fill_brewer(palette="Dark2") +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_blank(), #element_text(size=16, angle=20, hjust=1)
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=5),
    legend.title=element_text(size=5),
    panel.background = element_rect(fill = "white")
  ) +
  ggtitle("ACLAME - genes present at >10% abundance") +
  xlab('Sample ID') +
  scale_fill_tableau("Tableau 20", direction = -1) +
  ylab('Relative abundance') 

ggsave("Figure10_ICEberg_composition.jpeg", width = 30, height = 20, units = "cm")

ggplot(AMR_class_sum, aes(x = ID, y = sum_class, fill = Gene)) + 
  geom_bar(stat = "identity")+
  facet_wrap( ~ Study, scales='free_x',ncol = 2) +
  #scale_fill_brewer(palette="Dark2") +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_blank(), #element_text(size=16, angle=20, hjust=1)
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=5),
    legend.title=element_text(size=5),
    panel.background = element_rect(fill = "white")
  ) +
  ggtitle("       ACLAME - genes present at >10% abundance") +
  xlab('Sample ID') +
  scale_fill_tableau("Tableau 20", direction = -1) +
  ylab('CSS counts') 






##
# Figure 7 
##


AMR_class_sum_raw <- amr_melted_raw_analytic[Level_ID=="Class", .(sum_class= sum(Normalized_Count)),by=.(ID, Name, Study)]
AMR_class_sum_raw[,total:= sum(sum_class), by=.(ID)]
AMR_class_sum_raw[,percentage:= sum_class/total ,by=.(ID, Name) ]
AMR_class_sum_raw$Class <- AMR_class_sum_raw$Name
#AMR_class_sum[,percentage:= round(sum_class/total, digits=2) ,by=.(ID, Name) ] removes some with low proportions
ggplot(AMR_class_sum_raw, aes(x = ID, y = sum_class, fill = Class)) + 
  geom_bar(stat = "identity")+
  facet_wrap( ~ Study, scales='free_x',ncol = 2) +
  #scale_fill_brewer(palette="Dark2") +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_blank(), #element_text(size=16, angle=20, hjust=1)
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=10),
    legend.title=element_text(size=20),
    panel.background = element_rect(fill = "white")
  ) +
  xlab('Sample ID') +
  ylab('Raw mapped reads') 
