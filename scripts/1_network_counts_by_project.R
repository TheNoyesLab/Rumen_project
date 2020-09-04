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
Shi_samples_old_micro_names <- row.names(Shi_samples_phylum_microbiome_counts) 
#row.names(Shi_samples_phylum_microbiome_counts) <- paste("Shi_samples",Shi_samples_old_micro_names, sep=".")

Shi_samples_phylum_microbiome_counts[is.na(Shi_samples_phylum_microbiome_counts)] <- 0
t_Shi_samples_phylum_microbiome_counts <- t(Shi_samples_phylum_microbiome_counts)

Shi_Microbiome_Phylum.dt <- as.data.table(t_Shi_samples_phylum_microbiome_counts, keep.rownames = TRUE)
setkey(Shi_Microbiome_Phylum.dt, rn)

final_Shi_Microbiome_Phylum_counts <- analysis_metadata[Shi_Microbiome_Phylum.dt]
write.csv(final_Shi_Microbiome_Phylum_counts, "network_counts/counts_Shi_samples_phylummicro_aclame_iceberg_megaresclass.csv", row.names = FALSE)

#
## Dichotemize counts
#
dichotemized_Shi_Microbiome_Phylum.dt <- as.data.table(Shi_Microbiome_Phylum.dt)
dichotemized_Shi_Microbiome_Phylum.dt[, names(dichotemized_Shi_Microbiome_Phylum.dt)[-1] := lapply(.SD, function(x) as.integer(x!=0)), .SDcols = 2:ncol(dichotemized_Shi_Microbiome_Phylum.dt)]
setkey(dichotemized_Shi_Microbiome_Phylum.dt, rn)

final_Shi_Microbiome_Phylum_dichotemized_counts <- analysis_metadata[dichotemized_Shi_Microbiome_Phylum.dt]

write.csv(final_Shi_Microbiome_Phylum_dichotemized_counts, "network_counts/dicho_Shi_samples_phylummicro_aclame_iceberg_megaresclass.csv", row.names = FALSE)

###
#####
# Hess_samples
#####
###
Hess_samples_phylum_microbiome_samples = which(pData(microbiome_analytic_data[[2]])$Study == "Hess" )
Hess_samples_phylum_microbiome <- microbiome_analytic_data[[2]][, Hess_samples_phylum_microbiome_samples]
Hess_samples_phylum_microbiome_counts <- MRcounts(Hess_samples_phylum_microbiome)
Hess_samples_old_micro_names <- row.names(Hess_samples_phylum_microbiome_counts) 
#row.names(Hess_samples_phylum_microbiome_counts) <- paste("Hess_samples",Hess_samples_old_micro_names, sep=".")

Hess_samples_phylum_microbiome_counts[is.na(Hess_samples_phylum_microbiome_counts)] <- 0
t_Hess_samples_phylum_microbiome_counts <- t(Hess_samples_phylum_microbiome_counts)

Hess_Microbiome_Phylum.dt <- as.data.table(t_Hess_samples_phylum_microbiome_counts, keep.rownames = TRUE)
setkey(Hess_Microbiome_Phylum.dt, rn)

final_Hess_Microbiome_Phylum_counts <- analysis_metadata[Hess_Microbiome_Phylum.dt]
write.csv(final_Hess_Microbiome_Phylum_counts, "network_counts/counts_Hess_samples_phylummicro_aclame_iceberg_megaresclass.csv", row.names = FALSE)

#
## Dichotemize counts
#
dichotemized_Hess_Microbiome_Phylum.dt <- as.data.table(Hess_Microbiome_Phylum.dt)
dichotemized_Hess_Microbiome_Phylum.dt[, names(dichotemized_Hess_Microbiome_Phylum.dt)[-1] := lapply(.SD, function(x) as.integer(x!=0)), .SDcols = 2:ncol(dichotemized_Hess_Microbiome_Phylum.dt)]
setkey(dichotemized_Hess_Microbiome_Phylum.dt, rn)

final_Hess_Microbiome_Phylum_dichotemized_counts <- analysis_metadata[dichotemized_Hess_Microbiome_Phylum.dt]

write.csv(final_Hess_Microbiome_Phylum_dichotemized_counts, "network_counts/dicho_Hess_samples_phylummicro_aclame_iceberg_megaresclass.csv", row.names = FALSE)


###
#####
# Stewart_samples
#####
###
Stewart_samples_phylum_microbiome_samples = which(pData(microbiome_analytic_data[[2]])$Study == "Stewart" )
Stewart_samples_phylum_microbiome <- microbiome_analytic_data[[2]][, Stewart_samples_phylum_microbiome_samples]
Stewart_samples_phylum_microbiome_counts <- MRcounts(Stewart_samples_phylum_microbiome)
Stewart_samples_old_micro_names <- row.names(Stewart_samples_phylum_microbiome_counts) 
#row.names(Stewart_samples_phylum_microbiome_counts) <- paste("Stewart_samples",Stewart_samples_old_micro_names, sep=".")

Stewart_samples_phylum_microbiome_counts[is.na(Stewart_samples_phylum_microbiome_counts)] <- 0
t_Stewart_samples_phylum_microbiome_counts <- t(Stewart_samples_phylum_microbiome_counts)

Stewart_Microbiome_Phylum.dt <- as.data.table(t_Stewart_samples_phylum_microbiome_counts, keep.rownames = TRUE)
setkey(Stewart_Microbiome_Phylum.dt, rn)

final_Stewart_Microbiome_Phylum_counts <- analysis_metadata[Stewart_Microbiome_Phylum.dt]
write.csv(final_Stewart_Microbiome_Phylum_counts, "network_counts/counts_Stewart_samples_phylummicro_aclame_iceberg_megaresclass.csv", row.names = FALSE)

#
## Dichotemize counts
#
dichotemized_Stewart_Microbiome_Phylum.dt <- as.data.table(Stewart_Microbiome_Phylum.dt)
dichotemized_Stewart_Microbiome_Phylum.dt[, names(dichotemized_Stewart_Microbiome_Phylum.dt)[-1] := lapply(.SD, function(x) as.integer(x!=0)), .SDcols = 2:ncol(dichotemized_Stewart_Microbiome_Phylum.dt)]
setkey(dichotemized_Stewart_Microbiome_Phylum.dt, rn)

final_Stewart_Microbiome_Phylum_dichotemized_counts <- analysis_metadata[dichotemized_Stewart_Microbiome_Phylum.dt]

write.csv(final_Stewart_Microbiome_Phylum_dichotemized_counts, "network_counts/dicho_Stewart_samples_phylummicro_aclame_iceberg_megaresclass.csv", row.names = FALSE)

###
#####
# Wallace_samples
#####
###
Wallace_samples_phylum_microbiome_samples = which(pData(microbiome_analytic_data[[2]])$Study == "Wallace" )
Wallace_samples_phylum_microbiome <- microbiome_analytic_data[[2]][, Wallace_samples_phylum_microbiome_samples]
Wallace_samples_phylum_microbiome_counts <- MRcounts(Wallace_samples_phylum_microbiome)
Wallace_samples_old_micro_names <- row.names(Wallace_samples_phylum_microbiome_counts) 
#row.names(Wallace_samples_phylum_microbiome_counts) <- paste("Wallace_samples",Wallace_samples_old_micro_names, sep=".")

Wallace_samples_phylum_microbiome_counts[is.na(Wallace_samples_phylum_microbiome_counts)] <- 0
t_Wallace_samples_phylum_microbiome_counts <- t(Wallace_samples_phylum_microbiome_counts)

Wallace_Microbiome_Phylum.dt <- as.data.table(t_Wallace_samples_phylum_microbiome_counts, keep.rownames = TRUE)
setkey(Wallace_Microbiome_Phylum.dt, rn)

final_Wallace_Microbiome_Phylum_counts <- analysis_metadata[Wallace_Microbiome_Phylum.dt]
write.csv(final_Wallace_Microbiome_Phylum_counts, "network_counts/counts_Wallace_samples_phylummicro_aclame_iceberg_megaresclass.csv", row.names = FALSE)

#
## Dichotemize counts
#
dichotemized_Wallace_Microbiome_Phylum.dt <- as.data.table(Wallace_Microbiome_Phylum.dt)
dichotemized_Wallace_Microbiome_Phylum.dt[, names(dichotemized_Wallace_Microbiome_Phylum.dt)[-1] := lapply(.SD, function(x) as.integer(x!=0)), .SDcols = 2:ncol(dichotemized_Wallace_Microbiome_Phylum.dt)]
setkey(dichotemized_Wallace_Microbiome_Phylum.dt, rn)

final_Wallace_Microbiome_Phylum_dichotemized_counts <- analysis_metadata[dichotemized_Wallace_Microbiome_Phylum.dt]

write.csv(final_Wallace_Microbiome_Phylum_dichotemized_counts, "network_counts/dicho_Wallace_samples_phylummicro_aclame_iceberg_megaresclass.csv", row.names = FALSE)


