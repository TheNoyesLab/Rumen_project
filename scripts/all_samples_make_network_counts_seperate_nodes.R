# shi_samples
shi_samples_phylum_microbiome_samples = which(pData(microbiome_analytic_data[[2]])$Study == "Shi" )
shi_samples_phylum_microbiome <- microbiome_analytic_data[[2]][, shi_samples_phylum_microbiome_samples]
shi_samples_phylum_microbiome_counts <- MRcounts(shi_samples_phylum_microbiome)
shi_samples_old_micro_names <- row.names(shi_samples_phylum_microbiome_counts) 
row.names(shi_samples_phylum_microbiome_counts) <- paste("shi_samples",shi_samples_old_micro_names, sep=".")

# hess_samples
hess_samples_phylum_microbiome_samples = which(pData(microbiome_analytic_data[[2]])$Study == "Hess" )
hess_samples_phylum_microbiome <- microbiome_analytic_data[[2]][, hess_samples_phylum_microbiome_samples]
hess_samples_phylum_microbiome_counts <- MRcounts(hess_samples_phylum_microbiome)
hess_samples_old_micro_names <- row.names(hess_samples_phylum_microbiome_counts) 
row.names(hess_samples_phylum_microbiome_counts) <- paste("hess_samples",hess_samples_old_micro_names, sep=".")

# Wallace_samples
wallace_samples_phylum_microbiome_samples = which(pData(microbiome_analytic_data[[2]])$Study == "Wallace" )
wallace_samples_phylum_microbiome <- microbiome_analytic_data[[2]][, wallace_samples_phylum_microbiome_samples]
wallace_samples_phylum_microbiome_counts <- MRcounts(wallace_samples_phylum_microbiome)
wallace_samples_old_micro_names <- row.names(wallace_samples_phylum_microbiome_counts) 
row.names(wallace_samples_phylum_microbiome_counts) <- paste("wallace_samples",wallace_samples_old_micro_names, sep=".")

# stewart_samples
stewart_samples_phylum_microbiome_samples = which(pData(microbiome_analytic_data[[2]])$Study == "Stewart" )
stewart_samples_phylum_microbiome <- microbiome_analytic_data[[2]][, stewart_samples_phylum_microbiome_samples]
stewart_samples_phylum_microbiome_counts <- MRcounts(stewart_samples_phylum_microbiome)
stewart_samples_old_micro_names <- row.names(stewart_samples_phylum_microbiome_counts) 
row.names(stewart_samples_phylum_microbiome_counts) <- paste("stewart_samples",stewart_samples_old_micro_names, sep=".")



# Merge kuner samples
kuner_microbiome_phylum_feces_LA_labeled <- merge(kuner_LA_phylum_microbiome_counts, kuner_feces_phylum_microbiome_counts ,by = "row.names", all = TRUE)
kuner_microbiome_phylum_soil_labeled <- merge(kuner_soilsurface_phylum_microbiome_counts,kuner_soilcores_phylum_microbiome_counts ,by = "row.names", all = TRUE)
all_kuner_microbiome_phylum_labeled <- merge(kuner_microbiome_phylum_feces_LA_labeled,kuner_microbiome_phylum_soil_labeled ,by = "Row.names", all = TRUE)

# Merge XIT samples
XIT_microbiome_phylum_feces_LA_labeled <- merge(XIT_LA_phylum_microbiome_counts, XIT_feces_phylum_microbiome_counts ,by = "row.names", all = TRUE)
XIT_soilsurface_phylum_microbiome_counts <- dplyr::as_data_frame(XIT_soilsurface_phylum_microbiome_counts, rownames = "Row.names")
#XIT_soilsurface_phylum_microbiome_counts$Row.names <- row.names(XIT_soilsurface_phylum_microbiome_counts)
all_XIT_microbiome_phylum_feces_LA_labeled <- merge(XIT_soilsurface_phylum_microbiome_counts,XIT_microbiome_phylum_feces_LA_labeled, by = "Row.names", all = TRUE)

# Merge tylan 1
tylan1_phylum_feces_LA_labeled <- merge(tylan1_LA_phylum_microbiome_counts, tylan1_feces_phylum_microbiome_counts,by = "row.names", all = TRUE)
tylan1_soil_phylum_microbiome_counts <- dplyr::as_data_frame(tylan1_soil_phylum_microbiome_counts, rownames = "Row.names")
all_tylan1_microbiome_phylum_labeled <- merge(tylan1_soil_phylum_microbiome_counts,tylan1_phylum_feces_LA_labeled)

# combine all those counts
microbiome_phylum_labeled <- merge(all_XIT_microbiome_phylum_feces_LA_labeled,all_kuner_microbiome_phylum_labeled ,by = "Row.names", all = TRUE)
microbiome_phylum_labeled <- merge(microbiome_phylum_labeled, all_tylan1_microbiome_phylum_labeled,by = "Row.names", all = TRUE)

# Change this line if you only want a subset
#microbiome_phylum_labeled <- kuner_LA_phylum_microbiome_counts

microbiome_phylum_labeled[is.na(microbiome_phylum_labeled)] <- 0

# If only using one MRexperiment, don't need the following two commands
row.names(microbiome_phylum_labeled) <- microbiome_phylum_labeled$Row.names
microbiome_phylum_labeled <- microbiome_phylum_labeled[-1]



# Transpose counts
#AMR_Class <- t(AMR_Class)
Microbiome_Phylum <- t(microbiome_phylum_labeled)

# combine AMR and Microbiome counts
#combined_AMRMechanism_Phylum_counts <- merge(AMR_Class, Microbiome_Phylum, by="row.names", all = TRUE)
#set empty na columns to 0
Microbiome_Phylum[is.na(Microbiome_Phylum)] <- 0

#combined_AMRMechanism_Phylum_counts <- as.data.table(combined_AMRMechanism_Phylum_counts)

#setkey(combined_AMRMechanism_Phylum_counts, Row.names)
Microbiome_Phylum.dt <- as.data.table(Microbiome_Phylum, keep.rownames = TRUE)
setkey(Microbiome_Phylum.dt, rn)

#
## Dichotemize counts
#
dichotemized_combined_AMRMechanism_Phylum_counts <- Microbiome_Phylum.dt
dichotemized_combined_AMRMechanism_Phylum_counts[, names(dichotemized_combined_AMRMechanism_Phylum_counts)[-1] := lapply(.SD, function(x) as.integer(x!=0)), .SDcols = 2:ncol(dichotemized_combined_AMRMechanism_Phylum_counts)]
setkey(dichotemized_combined_AMRMechanism_Phylum_counts, rn)


# Make analysis metadata
# analysis_metadata <- microbiome_metadata[,.(ID, Group,Head,PREVCAT_A_APLUS, PREVCAT_ALL,microbiome_phylum_Richness,
#                                             microbiome_phylum_Shannon,microbiome_order_Richness,microbiome_order_Shannon)]  

analysis_metadata <- microbiome_metadata[,.(ID, Group,PREVCAT_A_APLUS, PREVCAT_ALL)]
setkey(analysis_metadata,ID)

#combined_Phylum_counts <- analysis_metadata[Microbiome_Phylum.dt]
combined_Phylum_counts <- analysis_metadata[dichotemized_combined_AMRMechanism_Phylum_counts]


#set empty na columns to 0
#combined_AMRMechanism_Phylum_counts[is.na(combined_AMRMechanism_Phylum_counts)] <- 0

#proj3_combined_AMRMechanism_Phylum_counts <- combined_AMRMechanism_Phylum_counts[Sample_type == "Individual"]
#proj4_combined_AMRMechanism_Phylum_counts <- combined_AMRMechanism_Phylum_counts[Sample_type == "Pen"]


#write to csv file, row.names=FALSE removes the 1st column
#write.csv(combined_Phylum_counts, "kuner_LA_Phylum_nodes.csv", row.names = FALSE)
write.csv(combined_Phylum_counts, "dicho_combined_tylan1_XIT_kuner_samples_Phylum_counts_seperate_nodes.csv", row.names = FALSE)

