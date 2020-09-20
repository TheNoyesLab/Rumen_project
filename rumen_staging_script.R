## Start with this staging file to set up your analysis.
# Source the utility functions file, which should be in the scripts folder with this file
source('scripts/meg_utility_functions.R')
source('scripts/load_libraries.R')

# Set working directory to the MEG_R_metagenomic_analysis folder and add your data to that folder
#setwd("")

# Set the output directory for graphs:
graph_output_dir = 'graphs'
# Set the output directory for statistics:
stats_output_dir = 'stats'
# In which column of the metadata file are the sample IDs stored?
sample_column_id = 'ID'

####################
## File locations ##
####################
## The files you want to use for input to this (for the MEG group analyses)
## is the AMR_analytic_matrix.csv. So you should have pulled these files from the output of the nextflow pipeline
## and you are now performing this analysis on your local machine.

## For the AMR analysis, you will also need to download the megares_annotations.csv
## file from the MEGARes website; the annotation file must be from the same version
## of the database as the file you used in the AmrPlusPlus pipeline, i.e. the headers
## must match between the annotation file and the database file.

amr_metadata_filepath = 'rumen_sample_metadata.csv'

###################
## User Controls ##
###################
## Hopefully, this section should be the only code you need to modify.
## However, you can look into the code in further sections if you need
## to change other, more subtle variables in the exploratory or
## statistical functions.

# The following is a list of analyses based on variables in
# your metadata.csv file that you want
# to use for EXPLORATORY analysis (NMDS, PCA, alpha rarefaction, barplots)
# NOTE: Exploratory variables cannot be numeric.


AMR_exploratory_analyses = list(
  # Analysis
  # Description:
  list(
    name = 'Study',
    subsets = list(),
    exploratory_var = 'Study',
    order = ''
  ),
  # Analysis
  # Description:
  list(
    name = 'Description',
    subsets = list(),
    exploratory_var = 'Description',
    order = ''
  ),
  # Analysis
  # Description:
  list(
    name = 'Species',
    subsets = list(),
    exploratory_var = 'Species',
    order = ''
  ),
  # Analysis
  # Description:
  list(
    name = 'Sequencing',
    subsets = list(),
    exploratory_var = 'Sequencing',
    order = ''
  ),
  # Shi Analysis
  # Description:
  list(
    name = 'Shi_Description',
    subsets = list('Study == Shi'),
    exploratory_var = 'Description',
    order = ''
  ),
  # Hess Analysis
  # Description:
  list(
    name = 'Hess_Description',
    subsets = list('Study == Hess'),
    exploratory_var = 'Description',
    order = ''
  ),
  # Stewart Analysis
  # Description:
  list(
    name = 'Stewart_Description',
    subsets = list('Study == Stewart'),
    exploratory_var = 'Description',
    order = ''
  ),
  # Wallace Analysis
  # Description:
  list(
    name = 'Wallace_Description',
    subsets = list('Study == Wallace'),
    exploratory_var = 'Description',
    order = ''
  )
)



microbiome_exploratory_analyses = list(
  # Analysis
  # Description:
  list(
    name = 'Study',
    subsets = list(),
    exploratory_var = 'Study',
    order = ''
  ),
  # Analysis
  # Description:
  list(
    name = 'Description',
    subsets = list(),
    exploratory_var = 'Description',
    order = ''
  ),
  # Analysis
  # Description:
  list(
    name = 'Species',
    subsets = list(),
    exploratory_var = 'Species',
    order = ''
  ),
  # Analysis
  # Description:
  list(
    name = 'Sequencing',
    subsets = list(),
    exploratory_var = 'Sequencing',
    order = ''
  ),
  # Shi Analysis
  # Description:
  list(
    name = 'Shi_Description',
    subsets = list('Study == Shi'),
    exploratory_var = 'Description',
    order = ''
  ),
  # Hess Analysis
  # Description:
  list(
    name = 'Hess_Description',
    subsets = list('Study == Hess'),
    exploratory_var = 'Description',
    order = ''
  ),
  # Stewart Analysis
  # Description:
  list(
    name = 'Stewart_Description',
    subsets = list('Study == Stewart'),
    exploratory_var = 'Description',
    order = ''
  ),
  # Wallace Analysis
  # Description:
  list(
    name = 'Wallace_Description',
    subsets = list('Study == Wallace'),
    exploratory_var = 'Description',
    order = ''
  )
)

# Each analyses you wish to perform should have its own list in the following
# statistical_analyses list.  A template is provided to get you started.
# Multiple analyses, subsets, and contrasts are valid, but only one random
# effect can be used per analysis.  The contrasts of interest must have their
# parent variable in the model matrix equation.  Contrasts are named by
# parent variable then child variable without a space inbetween, for example:
# PVar1Cvar1 where the model matrix equation is ~ 0 + Pvar1.

# This is not updated for the rumen data
AMR_statistical_analyses = list(
  # Analysis 1
  # Description:
  list(
    name = 'Treatment',
    subsets = list(),
    model_matrix = '~ 0 + Treatment ',
    contrasts = list('TreatmentCONV - TreatmentRWA'),
    random_effect = NA
  )
)

microbiome_statistical_analyses = list(
  # Analysis 1
  # Description:
  list(
    name = 'Treatment',
    subsets = list(),
    model_matrix = '~ 0 + Treatment ',
    contrasts = list('TreatmentCONV - TreatmentRWA'),
    random_effect = NA
  )
)

## Loading count tables

# Where is the metadata file stored on your machine?

# MEGARes
amr_count_matrix_filepath = 'MEGARes_AMR_analytic_matrix.csv'
# Name of the megares annotation file used for this project
megares_annotation_filename = 'data/amr/megares_modified_annotations_v2.0.csv'
source('scripts/metagenomeSeq_megaresv2.R')
# Save MRExperiments to another object
MEGARes_MRexperiment <- amr

# Running the following command creates exploratory figures and outputs count matrices
source('scripts/print_AMR_figures.R')
# NOTE!
##
### After creating the exploratory figures, you'll need to rename the following folders; "graphs/" and "amr_matrices"

# ACLAME (remove the step to remove sparse features)
amr_count_matrix_filepath = 'ACLAME_rumen_analytic_matrix.csv'
# Name of the megares annotation file used for this project
megares_annotation_filename = 'data/amr/aclame_db_annotations.csv'
source('scripts/metagenomeSeq_megaresv2.R')
ACLAME_MRexperiment <- amr
# NOTE!
##
### After creating the exploratory figures, you'll need to rename the following folders; "graphs/" and "amr_matrices"


## ICEBerg
amr_count_matrix_filepath = 'ICEBerg_AMR_analytic_matrix.csv'
# Name of the megares annotation file used for this project
megares_annotation_filename = 'data/amr/ICEBerg_db_annotations.csv'
source('scripts/metagenomeSeq_megaresv2.R')
ICEBerg_MRexperiment <- amr
# NOTE!
##
### After creating the exploratory figures, you'll need to rename the following folders; "graphs/" and "amr_matrices"

# Combined data, # this still causes an issue because "mergeMRexperiments()" does not merge rows by sample names
# and will instead make new labels for each sample
#allAMR_data <- mergeMRexperiments(MEGARes_MRexperiment, ACLAME_MRexperiment)
#allAMR_data <- mergeMRexperiments(allAMR_data, ICEBerg_MRexperiment)

#source('scripts/metagenomeSeq_combined_resistome.R')

#################################
## Microbiome - 16S or kraken? ##
#################################

# Where is the metadata file for the microbiome samples stored on your machine?
microbiome_temp_metadata_file = "rumen_sample_metadata.csv"

# If you used the AMR++ pipeline and have the kraken2 count matrix, point to the kraken file or see below for using qiime2 results.
kraken_temp_file = "kraken_analytic_matrix.csv"

## Run the analysis for the microbiome data
source('scripts/metagenomeSeq_kraken.R')

# After running this script, these are the useful objects that contain all the data aggregated to different levels
# The metagenomeSeq objects are contained in these lists "AMR_analytic_data" and "microbiome_analytic_data"
# Melted counts are contained in these data.table objects "amr_melted_analytic" "microbiome_melted_analytic"

## Run code to make some exploratory figures and output count matrices.
# These are commented out for now because they take a long time to run.
source('scripts/print_microbiome_figures.R')



### Start of code for extra figures

setkey(amr_melted_raw_analytic,ID)
setkey(amr_melted_analytic,ID)

setkey(microbiome_melted_analytic,ID)
# Set keys for both metadata files
setkey(metadata,ID)
setkey(microbiome_metadata,ID)
microbiome_melted_analytic <- microbiome_melted_analytic[microbiome_metadata]
amr_melted_raw_analytic <- amr_melted_raw_analytic[metadata]
amr_melted_analytic <- amr_melted_analytic[metadata]


## Microbiome counts
microbiome_phylum_sum <- microbiome_melted_analytic[Level_ID == "Phylum", .(sum_phylum= sum(Normalized_Count)),by=.(ID, Name, Level_ID,Study)]
microbiome_phylum_sum[,total:= sum(sum_phylum), by=.(ID)]
microbiome_phylum_sum[,percentage:= sum_phylum/total ,by=.(ID, Name) ]
# only keep taxa greater than 1%
microbiome_phylum_sum <- microbiome_phylum_sum[percentage > .01]
microbiome_phylum_sum[,total:= sum(sum_phylum), by=.(ID)]
microbiome_phylum_sum[,percentage:= sum_phylum/total ,by=.(ID, Name) ]

microbiome_phylum_sum$Name = droplevels(microbiome_phylum_sum$Name)
#microbiome_phylum_sum$Name = factor(microbiome_phylum_sum$Name ,levels=c("Planctomycetes","Gemmatimonadetes","Chloroflexi","Verrucomicrobia", "Acidobacteria","Actinobacteria" ,"Bacteroidetes" , "Proteobacteria" , "Firmicutes"))


microbiome_phylum_sum$Phylum <- microbiome_phylum_sum$Name
#microbiome_phylum_sum[,percentage:= round(sum_class/total, digits=2) ,by=.(ID, Name) ] removes some with low proportions

ggplot(microbiome_phylum_sum, aes(x = ID, y = percentage, fill = Phylum)) +
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
    legend.text=element_text(size=10),
    legend.title=element_text(size=20),
    panel.background = element_rect(fill = "white")
  ) +
  scale_fill_tableau("Tableau 20") +
  ggtitle("Microbiome composition in by treatment (only taxa > 1% per sample)") +
  xlab('Sample ID') +
  ylab('Relative abundance')
