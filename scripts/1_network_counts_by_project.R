# metadata
analysis_metadata <- microbiome_metadata[,.(ID, Description, Species)]
setkey(analysis_metadata,ID)

###
#####
# Shi_samples
#####
###
Shi_samples_phylum_microbiome_samples = which(pData(microbiome_analytic_data[[2]])$Study == "Shi" )
Shi_samples_phylum_microbiome <- microbiome_analytic_data[[2]][, Shi_samples_phylum_microbiome_samples]
Shi_samples_phylum_microbiome_counts <- MRcounts(Shi_samples_phylum_microbiome)
Shi_samples_old_microbiome_names <- row.names(Shi_samples_phylum_microbiome_counts) 
row.names(Shi_samples_phylum_microbiome_counts) <- paste("micro",Shi_samples_old_microbiome_names, sep=".")

Shi_samples_ACLAME_samples = which(pData(ACLAME_analytic_data[[1]])$Study == "Shi" )
Shi_samples_ACLAME <- ACLAME_analytic_data[[1]][, Shi_samples_ACLAME_samples]
Shi_samples_ACLAME_counts <- MRcounts(Shi_samples_ACLAME)
Shi_samples_old_ACLAME_names <- row.names(Shi_samples_ACLAME_counts) 
row.names(Shi_samples_ACLAME_counts) <- paste("ACLAME",Shi_samples_old_ACLAME_names, sep=".")

Shi_samples_MEGARes_samples = which(pData(MEGARes_analytic_data[[1]])$Study == "Shi" )
Shi_samples_MEGARes <- MEGARes_analytic_data[[1]][, Shi_samples_MEGARes_samples]
Shi_samples_MEGARes_counts <- MRcounts(Shi_samples_MEGARes)
Shi_samples_old_MEGARes_names <- row.names(Shi_samples_MEGARes_counts) 
row.names(Shi_samples_MEGARes_counts) <- paste("AMR",Shi_samples_old_MEGARes_names, sep=".")

Shi_samples_ICEBerg_samples = which(pData(ICEBerg_analytic_data[[1]])$Study == "Shi" )
Shi_samples_ICEBerg <- ICEBerg_analytic_data[[1]][, Shi_samples_ICEBerg_samples]
Shi_samples_ICEBerg_counts <- MRcounts(Shi_samples_ICEBerg)
Shi_samples_old_ICEBerg_names <- row.names(Shi_samples_ICEBerg_counts) 
row.names(Shi_samples_ICEBerg_counts) <- paste("ICE",Shi_samples_old_ICEBerg_names, sep=".")


# Combine counts ( need to be transposed)
Shi_micro_aclame_counts  <- merge(t(Shi_samples_phylum_microbiome_counts), t(Shi_samples_ACLAME_counts),by = "row.names", all.x = TRUE) 
Shi_megares_ICE_counts <-  merge(t(Shi_samples_MEGARes_counts), t(Shi_samples_ICEBerg_counts),by = "row.names", all.x = TRUE) 

Shi_all_counts <- merge(Shi_megares_ICE_counts, Shi_micro_aclame_counts, by = "Row.names", all = TRUE)

# Replace NA with 0 and transpose
Shi_all_counts[is.na(Shi_all_counts)] <- 0

# make dt, add metadata
Shi_all_counts.dt <- as.data.table(Shi_all_counts, keep.rownames = FALSE)
setkey(Shi_all_counts.dt, Row.names)
final_Shi_network_counts <- analysis_metadata[Shi_all_counts.dt]

# write the counts
write.csv(final_Shi_network_counts, "network_counts/counts_Shi_samples_phylummicro_aclame_iceberg_megaresclass.csv", row.names = FALSE)

#
## Dichotemize counts
#
dichotemized_Shi_all_counts.dt <- as.data.table(Shi_all_counts.dt)
dichotemized_Shi_all_counts.dt[, names(dichotemized_Shi_all_counts.dt)[-1] := lapply(.SD, function(x) as.integer(x!=0)), .SDcols = 2:ncol(dichotemized_Shi_all_counts.dt)]
setkey(dichotemized_Shi_all_counts.dt, Row.names)

final_Shi_dichotemized_network_counts <- analysis_metadata[dichotemized_Shi_all_counts.dt]

write.csv(final_Shi_dichotemized_network_counts, "network_counts/dicho_Shi_samples_phylummicro_aclame_iceberg_megaresclass.csv", row.names = FALSE)

###
#####
# Hess_samples
#####
###
Hess_samples_phylum_microbiome_samples = which(pData(microbiome_analytic_data[[2]])$Study == "Hess" )
Hess_samples_phylum_microbiome <- microbiome_analytic_data[[2]][, Hess_samples_phylum_microbiome_samples]
Hess_samples_phylum_microbiome_counts <- MRcounts(Hess_samples_phylum_microbiome)
Hess_samples_old_microbiome_names <- row.names(Hess_samples_phylum_microbiome_counts) 
row.names(Hess_samples_phylum_microbiome_counts) <- paste("micro",Hess_samples_old_microbiome_names, sep=".")

Hess_samples_ACLAME_samples = which(pData(ACLAME_analytic_data[[1]])$Study == "Hess" )
Hess_samples_ACLAME <- ACLAME_analytic_data[[1]][, Hess_samples_ACLAME_samples]
Hess_samples_ACLAME_counts <- MRcounts(Hess_samples_ACLAME)
Hess_samples_old_ACLAME_names <- row.names(Hess_samples_ACLAME_counts) 
row.names(Hess_samples_ACLAME_counts) <- paste("ACLAME",Hess_samples_old_ACLAME_names, sep=".")

Hess_samples_MEGARes_samples = which(pData(MEGARes_analytic_data[[1]])$Study == "Hess" )
Hess_samples_MEGARes <- MEGARes_analytic_data[[1]][, Hess_samples_MEGARes_samples]
Hess_samples_MEGARes_counts <- MRcounts(Hess_samples_MEGARes)
Hess_samples_old_MEGARes_names <- row.names(Hess_samples_MEGARes_counts) 
row.names(Hess_samples_MEGARes_counts) <- paste("AMR",Hess_samples_old_MEGARes_names, sep=".")

Hess_samples_ICEBerg_samples = which(pData(ICEBerg_analytic_data[[1]])$Study == "Hess" )
Hess_samples_ICEBerg <- ICEBerg_analytic_data[[1]][, Hess_samples_ICEBerg_samples]
Hess_samples_ICEBerg_counts <- MRcounts(Hess_samples_ICEBerg)
Hess_samples_old_ICEBerg_names <- row.names(Hess_samples_ICEBerg_counts) 
row.names(Hess_samples_ICEBerg_counts) <- paste("ICE",Hess_samples_old_ICEBerg_names, sep=".")


# Combine counts ( need to be transposed)
Hess_micro_aclame_counts  <- merge(t(Hess_samples_phylum_microbiome_counts), t(Hess_samples_ACLAME_counts),by = "row.names", all.x = TRUE) 
Hess_megares_ICE_counts <-  merge(t(Hess_samples_MEGARes_counts), t(Hess_samples_ICEBerg_counts),by = "row.names", all.x = TRUE) 

Hess_all_counts <- merge(Hess_megares_ICE_counts, Hess_micro_aclame_counts, by = "Row.names", all = TRUE)

# Replace NA with 0 and transpose
Hess_all_counts[is.na(Hess_all_counts)] <- 0

# make dt, add metadata
Hess_all_counts.dt <- as.data.table(Hess_all_counts, keep.rownames = FALSE)
setkey(Hess_all_counts.dt, Row.names)
final_Hess_network_counts <- analysis_metadata[Hess_all_counts.dt]

# write the counts
write.csv(final_Hess_network_counts, "network_counts/counts_Hess_samples_phylummicro_aclame_iceberg_megaresclass.csv", row.names = FALSE)

#
## Dichotemize counts
#
dichotemized_Hess_all_counts.dt <- as.data.table(Hess_all_counts.dt)
dichotemized_Hess_all_counts.dt[, names(dichotemized_Hess_all_counts.dt)[-1] := lapply(.SD, function(x) as.integer(x!=0)), .SDcols = 2:ncol(dichotemized_Hess_all_counts.dt)]
setkey(dichotemized_Hess_all_counts.dt, Row.names)

final_Hess_dichotemized_network_counts <- analysis_metadata[dichotemized_Hess_all_counts.dt]

write.csv(final_Hess_dichotemized_network_counts, "network_counts/dicho_Hess_samples_phylummicro_aclame_iceberg_megaresclass.csv", row.names = FALSE)

###
#####
# Stewart_samples
#####
###
Stewart_samples_phylum_microbiome_samples = which(pData(microbiome_analytic_data[[2]])$Study == "Stewart" )
Stewart_samples_phylum_microbiome <- microbiome_analytic_data[[2]][, Stewart_samples_phylum_microbiome_samples]
Stewart_samples_phylum_microbiome_counts <- MRcounts(Stewart_samples_phylum_microbiome)
Stewart_samples_old_microbiome_names <- row.names(Stewart_samples_phylum_microbiome_counts) 
row.names(Stewart_samples_phylum_microbiome_counts) <- paste("micro",Stewart_samples_old_microbiome_names, sep=".")

Stewart_samples_ACLAME_samples = which(pData(ACLAME_analytic_data[[1]])$Study == "Stewart" )
Stewart_samples_ACLAME <- ACLAME_analytic_data[[1]][, Stewart_samples_ACLAME_samples]
Stewart_samples_ACLAME_counts <- MRcounts(Stewart_samples_ACLAME)
Stewart_samples_old_ACLAME_names <- row.names(Stewart_samples_ACLAME_counts) 
row.names(Stewart_samples_ACLAME_counts) <- paste("ACLAME",Stewart_samples_old_ACLAME_names, sep=".")

Stewart_samples_MEGARes_samples = which(pData(MEGARes_analytic_data[[1]])$Study == "Stewart" )
Stewart_samples_MEGARes <- MEGARes_analytic_data[[1]][, Stewart_samples_MEGARes_samples]
Stewart_samples_MEGARes_counts <- MRcounts(Stewart_samples_MEGARes)
Stewart_samples_old_MEGARes_names <- row.names(Stewart_samples_MEGARes_counts) 
row.names(Stewart_samples_MEGARes_counts) <- paste("AMR",Stewart_samples_old_MEGARes_names, sep=".")

Stewart_samples_ICEBerg_samples = which(pData(ICEBerg_analytic_data[[1]])$Study == "Stewart" )
Stewart_samples_ICEBerg <- ICEBerg_analytic_data[[1]][, Stewart_samples_ICEBerg_samples]
Stewart_samples_ICEBerg_counts <- MRcounts(Stewart_samples_ICEBerg)
Stewart_samples_old_ICEBerg_names <- row.names(Stewart_samples_ICEBerg_counts) 
row.names(Stewart_samples_ICEBerg_counts) <- paste("ICE",Stewart_samples_old_ICEBerg_names, sep=".")


# Combine counts ( need to be transposed)
Stewart_micro_aclame_counts  <- merge(t(Stewart_samples_phylum_microbiome_counts), t(Stewart_samples_ACLAME_counts),by = "row.names", all.x = TRUE) 
Stewart_megares_ICE_counts <-  merge(t(Stewart_samples_MEGARes_counts), t(Stewart_samples_ICEBerg_counts),by = "row.names", all.x = TRUE) 

Stewart_all_counts <- merge(Stewart_megares_ICE_counts, Stewart_micro_aclame_counts, by = "Row.names", all = TRUE)

# Replace NA with 0 and transpose
Stewart_all_counts[is.na(Stewart_all_counts)] <- 0

# make dt, add metadata
Stewart_all_counts.dt <- as.data.table(Stewart_all_counts, keep.rownames = FALSE)
setkey(Stewart_all_counts.dt, Row.names)
final_Stewart_network_counts <- analysis_metadata[Stewart_all_counts.dt]

# write the counts
write.csv(final_Stewart_network_counts, "network_counts/counts_Stewart_samples_phylummicro_aclame_iceberg_megaresclass.csv", row.names = FALSE)

#
## Dichotemize counts
#
dichotemized_Stewart_all_counts.dt <- as.data.table(Stewart_all_counts.dt)
dichotemized_Stewart_all_counts.dt[, names(dichotemized_Stewart_all_counts.dt)[-1] := lapply(.SD, function(x) as.integer(x!=0)), .SDcols = 2:ncol(dichotemized_Stewart_all_counts.dt)]
setkey(dichotemized_Stewart_all_counts.dt, Row.names)

final_Stewart_dichotemized_network_counts <- analysis_metadata[dichotemized_Stewart_all_counts.dt]

write.csv(final_Stewart_dichotemized_network_counts, "network_counts/dicho_Stewart_samples_phylummicro_aclame_iceberg_megaresclass.csv", row.names = FALSE)

###
#####
# Wallace_samples
#####
###
Wallace_samples_phylum_microbiome_samples = which(pData(microbiome_analytic_data[[2]])$Study == "Wallace" )
Wallace_samples_phylum_microbiome <- microbiome_analytic_data[[2]][, Wallace_samples_phylum_microbiome_samples]
Wallace_samples_phylum_microbiome_counts <- MRcounts(Wallace_samples_phylum_microbiome)
Wallace_samples_old_microbiome_names <- row.names(Wallace_samples_phylum_microbiome_counts) 
row.names(Wallace_samples_phylum_microbiome_counts) <- paste("micro",Wallace_samples_old_microbiome_names, sep=".")

Wallace_samples_ACLAME_samples = which(pData(ACLAME_analytic_data[[1]])$Study == "Wallace" )
Wallace_samples_ACLAME <- ACLAME_analytic_data[[1]][, Wallace_samples_ACLAME_samples]
Wallace_samples_ACLAME_counts <- MRcounts(Wallace_samples_ACLAME)
Wallace_samples_old_ACLAME_names <- row.names(Wallace_samples_ACLAME_counts) 
row.names(Wallace_samples_ACLAME_counts) <- paste("ACLAME",Wallace_samples_old_ACLAME_names, sep=".")

Wallace_samples_MEGARes_samples = which(pData(MEGARes_analytic_data[[1]])$Study == "Wallace" )
Wallace_samples_MEGARes <- MEGARes_analytic_data[[1]][, Wallace_samples_MEGARes_samples]
Wallace_samples_MEGARes_counts <- MRcounts(Wallace_samples_MEGARes)
Wallace_samples_old_MEGARes_names <- row.names(Wallace_samples_MEGARes_counts) 
row.names(Wallace_samples_MEGARes_counts) <- paste("AMR",Wallace_samples_old_MEGARes_names, sep=".")

Wallace_samples_ICEBerg_samples = which(pData(ICEBerg_analytic_data[[1]])$Study == "Wallace" )
Wallace_samples_ICEBerg <- ICEBerg_analytic_data[[1]][, Wallace_samples_ICEBerg_samples]
Wallace_samples_ICEBerg_counts <- MRcounts(Wallace_samples_ICEBerg)
Wallace_samples_old_ICEBerg_names <- row.names(Wallace_samples_ICEBerg_counts) 
row.names(Wallace_samples_ICEBerg_counts) <- paste("ICE",Wallace_samples_old_ICEBerg_names, sep=".")


# Combine counts ( need to be transposed)
Wallace_micro_aclame_counts  <- merge(t(Wallace_samples_phylum_microbiome_counts), t(Wallace_samples_ACLAME_counts),by = "row.names", all.x = TRUE) 
Wallace_megares_ICE_counts <-  merge(t(Wallace_samples_MEGARes_counts), t(Wallace_samples_ICEBerg_counts),by = "row.names", all.x = TRUE) 

Wallace_all_counts <- merge(Wallace_megares_ICE_counts, Wallace_micro_aclame_counts, by = "Row.names", all = TRUE)

# Replace NA with 0 and transpose
Wallace_all_counts[is.na(Wallace_all_counts)] <- 0

# make dt, add metadata
Wallace_all_counts.dt <- as.data.table(Wallace_all_counts, keep.rownames = FALSE)
setkey(Wallace_all_counts.dt, Row.names)
final_Wallace_network_counts <- analysis_metadata[Wallace_all_counts.dt]

# write the counts
write.csv(final_Wallace_network_counts, "network_counts/counts_Wallace_samples_phylummicro_aclame_iceberg_megaresclass.csv", row.names = FALSE)

#
## Dichotemize counts
#
dichotemized_Wallace_all_counts.dt <- as.data.table(Wallace_all_counts.dt)
dichotemized_Wallace_all_counts.dt[, names(dichotemized_Wallace_all_counts.dt)[-1] := lapply(.SD, function(x) as.integer(x!=0)), .SDcols = 2:ncol(dichotemized_Wallace_all_counts.dt)]
setkey(dichotemized_Wallace_all_counts.dt, Row.names)

final_Wallace_dichotemized_network_counts <- analysis_metadata[dichotemized_Wallace_all_counts.dt]

write.csv(final_Wallace_dichotemized_network_counts, "network_counts/dicho_Wallace_samples_phylummicro_aclame_iceberg_megaresclass.csv", row.names = FALSE)
