###############
## Functions ##
###############
## Utility functions that can be optionally used later for producing
## graphs and performing reshaping operations
## With some edits for a higher trymax and change to heatmap code

# Call variables from parent scope
`..` <- function (..., .env = sys.parent(2)) {
  get(deparse(substitute(...)), env = .env)
}

# Misc reshape function for data table
melt_dt <- function(D, level_id) {
  temp <- reshape2::melt(D, variable.name='Sample', value.name='Normalized_Count')
  names(temp) <- c('Name', 'ID', 'Normalized_Count')
  temp <- data.table(cbind(rep(level_id, nrow(temp)), temp))
  names(temp)[1] <- 'Level_ID'
  return(temp)
}

# Filter data by quantile
meg_filter_data <- function(data_list,
                            filter_min_threshold) {
  local_obj <- data_list
  for( l in 1:length(local_obj) ) {
    filter_threshold <- quantile(rowSums(MRcounts(local_obj[[l]])), 0.15)
    if( filter_threshold > filter_min_threshold ) filter_threshold <- filter_min_threshold
    local_obj[[l]] <- local_obj[[l]][which(rowSums(MRcounts(local_obj[[l]])) >= filter_threshold ), ]
    cumNorm(local_obj[[l]])
  }
  return(local_obj)
}




make_sparse <- function(df, rownames, excludes, filter_min_threshold=0.15) {
  local_df <- df[, .SD, .SDcols=!excludes]
  filter_threshold <- quantile(rowSums(local_df), 0.15)
  if( filter_threshold > filter_min_threshold ) filter_threshold <- filter_min_threshold
  chosen <- which(rowSums(local_df) >= filter_threshold )
  ret <- as.data.frame(local_df[chosen, ])
  rownames(ret) <- df[[rownames]][chosen]
  return(ret)
}

data_subset <- function(data_obj, subsets) {
  local_meta <- data.table(pData(data_obj))
  local_subset <- c()
  for( c in 1:length(subsets) ) {
    
    conditional_terms <- unlist(strsplit(subsets[[c]], ' '))
    conditional_string <- paste('local_meta[[\'', conditional_terms[1],
                                '\']] ', conditional_terms[2],
                                ' \'', conditional_terms[3], '\'',
                                sep='', collapse='')
    if(length(local_subset) > 0) {
      local_subset <- intersect(local_subset, which(eval(parse(text=conditional_string))))
    }
    else {
      local_subset <- which(eval(parse(text=conditional_string)))
    }
  }
  return(data_obj[, local_subset])
}


data_subset_long <- function(data_obj, subsets) {
  local_subset <- c()
  for( c in 1:length(subsets) ) {
    
    conditional_terms <- unlist(strsplit(subsets[[c]], ' '))
    conditional_string <- paste('data_obj[[\'', conditional_terms[1],
                                '\']] ', conditional_terms[2],
                                ' \'', conditional_terms[3], '\'',
                                sep='', collapse='')
    if(length(local_subset) > 0) {
      local_subset <- intersect(local_subset, which(eval(parse(text=conditional_string))))
    }
    else {
      local_subset <- which(eval(parse(text=conditional_string)))
    }
  }
  return(data_obj[local_subset, ])
}


##############################
###                        ###
####     Ordination       ####
###                        ###
##############################
# Calculate the points for convex hulls around data for ordination plots
meg_find_hulls <- function(x) x[chull(x$Ord1, x$Ord2),]

hullPlot <- function(df, grouping, legend = TRUE){
  # requires that your columns are called NMDS1, NMDS2, NMDS3
  
  # requires that the "grouping" column is printed as per the name of the column
  # you want separate hulls for
  
  
  allSites <- sort(as.vector(unique(df[[grouping]])))
  matList <- list()
  hullList <- list()
  cols <- as.character(brewer_pal(type = "qual", palette = 1, direction = 1)(length(allSites)))
  
  
  # this loop creates the points for each site
  # it also calculates a separate hull for each site
  # "site" or whatever grouping variable you are using
  for(thisSite in 1:length(allSites)){
    
    tmp <- df[df[grouping] == allSites[thisSite], ]
    
    plot3d(tmp$NMDS1, tmp$NMDS2, tmp$NMDS3, col = cols[thisSite], box = FALSE,
           type = "s", radius = 0.01, add = ifelse(thisSite > 1, TRUE, FALSE),
           xlab = "", ylab = "", zlab = "")
    
    
    matList[[thisSite]] <- matrix(
      c(tmp[[grep("NMDS1", names(tmp))]],
        tmp[[grep("NMDS2", names(tmp))]],
        tmp[[grep("NMDS3", names(tmp))]]), ncol = 3)
    
    hullList[[thisSite]] <- t(convhulln(matList[[thisSite]]))
    
  }
  
  
  # this will run if you have legend = TRUE (the default)
  if(legend){
    # this slows down the plotting which is necessary otherwise
    # the printing can lag and the legend goes to a funny size
    Sys.sleep(0.5)
    
    # you can change your legend as per 'legend' commands
    legend3d("bottom", legend = allSites,
             # uses the same cols as for plotting
             col = cols,
             # symbol size
             pch = 16,
             inset=c(0.02),
             horiz = TRUE)
    
  }
  
  
  # this loop plots the hulls
  for(hull in seq_along(matList)){
    rgl.triangles(matList[[hull]][hullList[[hull]],1],matList[[hull]][hullList[[hull]],2],matList[[hull]][hullList[[hull]],3],
                  col=cols[hull],
                  
                  # change the alpha to change how see through they are
                  alpha=.6)
    
  }
  
}
# Function for computing ordination plots with convex hulls
# using ggplot2.  You will have to specify the facet variable,
# i.e. which experimental design variable determines the grouping of points.
#
# Function arguments:
#   data_list: a list containing the MRexperiment objects for all levels
#   data_names: a character vector of level names for each MRexperiment in data_list
#   metadata: a data.table of metadata for each sample
#   sample_var: the column name in metadata that specifies the sample IDs (must match MRexperiment columns)
#   hull_var: the metadata column name by which to group data points
#   outdir: the file path of the directory for output of files
#   data_type: the data type being computed on, e.g. AMR or Microbiome
#   method: choice of 'NMDS' or 'PCA' for method of ordination
#
#   outputs:
meg_ordination <- function(data_list,
                           data_names,
                           metadata,
                           sample_var,
                           hull_var,
                           analysis_subset,
                           outdir,
                           data_type,
                           method='NMDS',
                           factor_order) {
  all_ord <- data.table(ID=character(),
                        Level_ID=character(),
                        Ord1=numeric(),
                        Ord2=numeric(),
                        Group_Var=character())
  all_hulls <- data.table(Group_Var=character(),
                          Ord1=numeric(),
                          Ord2=numeric(),
                          Level_ID=character())
  setkey(all_ord, ID)
  
  local_obj <- data_list
  for( l in 1:length(local_obj) ) {
    
    if(length(analysis_subset) > 0) {
      local_obj[[l]] <- data_subset(local_obj[[l]], analysis_subset)
    }
    local_meta <- metadata
    local_meta <- local_meta[local_meta[[sample_var]] %in% colnames(MRcounts(local_obj[[l]])), ]
    
    # Transpose the matrix for NMDS (groups are now in rows and features in columns)
    t_data <- t(MRcounts(local_obj[[l]]))
    local_meta <- local_meta[rowSums(t_data) > 0, ]
    t_data <- t_data[rowSums(t_data) > 0, ]
    t_data <- t_data[which(!is.na(local_meta[[sample_var]]) & local_meta[[sample_var]] != 'NA'), colSums(t_data) > 0]
    
    if( method == 'NMDS' ) {
      # Set parallel to whatever your computer can support in terms of CPU count
      #ord.res <- metaMDS(t_data, autotransform=F, parallel=7, trymax=1000)  #### adding way more try max, from 49 to 1000
      ord.res <- metaMDS(t_data, autotransform=F, trymax=100, stratmax = 0.999999)  #### adding way more try max, from 49 to 1000, removed parallel=7,
      ord_points <- data.table(ord.res$points)
      names(ord_points) <- c('Ord1', 'Ord2')
      ord_points[, ID :=( rownames(ord.res$points) )]
    }
    else if( method == 'PCA' ) {
      ord.res <- prcomp(t_data, center=T, scale=T)
      
      # Format to include metadata for ggplot2
      ord_points <- data.table(ord.res$x[, 1:2])
      names(ord_points) <- c('Ord1', 'Ord2')
      ord_points[, ID :=( rownames(ord.res$x) )]
    }
    else {
      stop('method must be either NMDS or PCA')
    }
    
    ord_points[, Level_ID :=( rep(data_names[l], nrow(ord_points)) )]
    setkey(ord_points, ID)
    ord_points <- metadata[ord_points]
    ord_points <- ord_points[, .SD, .SDcols=c(sample_var, hull_var, 'Ord1', 'Ord2', 'Level_ID')]
    ## Check if needing to change factor order. Need to convert to data frame first, then back to data table
    ord_points <- as.data.frame(ord_points)
    ifelse(factor_order != '', ord_points[[hull_var]] <- factor(ord_points[[hull_var]], levels = eval(factor_order)),1)
    ord_points <- as.data.table(ord_points)
    names(ord_points)[2] <- 'Group_Var'
    all_ord <- rbind(all_ord, ord_points)
    
    hulls <- ord_points[, meg_find_hulls(.SD), .SDcols=c('Ord1', 'Ord2'), by=Group_Var]
    hulls[, Level_ID :=( rep(data_names[l], nrow(hulls)) )]
    hulls[, .SD, .SDcols=!'Level_ID']
    all_hulls <- rbind(all_hulls, hulls)
    
    # Plot graphs with convex hulls
    g_ord <- ggplot(data=ord_points, aes(Ord1, Ord2, color=Group_Var, fill=Group_Var)) +
      geom_point(size=2.5) + geom_polygon(data=hulls, aes(x=Ord1, y=Ord2, color=Group_Var, fill=Group_Var),
                                          alpha=0.2, show.legend=F)
    #g_ord <- g_ord + scale_fill_tableau("Classic 20", direction = -1)
    #g_ord <- g_ord + scale_color_tableau("Classic 20", direction = -1)
    g_ord <- g_ord +
      ggtitle(paste(method, ' for ', data_type, ' ', 'by ', hull_var, '\nAnnotation Level: ',
                    data_names[l],
                    sep='',
                    collapse='')) +
      labs(color=hull_var) +
      guides(fill=F) +
      theme(strip.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.text.x=element_blank(),
            axis.title.x=element_text(size=24),
            axis.title.y=element_text(size=24, hjust=0.5),
            #legend.position="right",
            legend.title=element_text(size=20),
            legend.text=element_text(size=18),
            plot.title=element_text(size=24, hjust=0.5),
            panel.background = element_rect(fill = "white"),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.spacing=unit(0.1, "lines"),
            panel.border = element_rect(colour = "black", fill=NA, size=0.5))
    if( method == 'NMDS' ) {
      g_ord <- g_ord + xlab('MDS1') + ylab('MDS2')
    }
    else if( method == 'PCA' ) {
      g_ord <- g_ord + xlab('PC1') + ylab('PC2')
    }
    
    # Open the graphics device at the specified location and figure size
    png(filename=paste(outdir, '/', method, '_', hull_var, '_',
                       data_names[l], '.png', sep='', collapse=''),
        width=1024, height=768)
    
    print(g_ord)
    
    # Turn off graphics device to save the graphic
    dev.off()
    
  }
  
  
  all_ord <- within(all_ord, Level_ID
                    <- factor(Level_ID, levels=data_names,
                              ordered=T))
  all_hulls <- within(all_hulls, Level_ID
                      <- factor(Level_ID, levels=data_names,
                                ordered=T))
  
  g_all_ord <- ggplot(data=all_ord, aes(Ord1, Ord2, color=Group_Var, fill=Group_Var)) +
    geom_point(size=3) +
    geom_polygon(data=all_hulls, aes(x=Ord1, y=Ord2, color=Group_Var, fill=Group_Var),
                 alpha=0.2, show.legend=F) +
    facet_wrap(~Level_ID, nrow=2)
  g_all_ord <- g_all_ord +
    ggtitle(paste(method, ' for ', data_type, ' by ', hull_var, sep='', collapse='')) +
    labs(color=hull_var) +
    guides(fill=F) +
    theme(strip.text.x=element_text(size=26),
          axis.text.y=element_blank(),
          axis.text.x=element_blank(),
          axis.title.x=element_text(size=26),
          axis.title.y=element_text(size=26, hjust=0.5),
          #legend.position="right",
          legend.title=element_text(size=24, hjust=0.5),
          legend.text=element_text(size=20),
          plot.title=element_text(size=30, hjust=0.5),
          panel.background = element_rect(fill = "white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing=unit(0, "mm"),
          #strip.background = element_rect(size = 2, colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))
  #g_all_ord <- g_all_ord + scale_color_tableau("Classic 20", direction = -1)
  #g_all_ord <- g_all_ord + scale_fill_tableau("Classic 20", direction = -1)
  
  png(filename=paste(outdir, '/', method, '_', hull_var, '_',
                     'AllLevels.png', sep='', collapse=''),
      width=1024, height=768)
  print(g_all_ord)
  dev.off()
  return(all_ord)
}

##############################
###                        ###
####      Heatmap         ####
###                        ###
##############################
heatmap_select_top_counts <- function(X, group_var, sample_var, n) {
  return(X[, tail(.SD, n), by=c(group_var, sample_var)])
}

heatmap_col <- c("#849bed", "#9ce2de","#f7f6ed","#efdfa7","#f77a74")

#heatmap_col <- c("grey","#f7f6ed","#efdfa7","#f77a74")

meg_heatmap <- function(melted_data,
                        metadata,
                        sample_var,
                        group_var,
                        level_var,
                        analysis_subset,
                        outdir,
                        data_type,
                        factor_order='') {
  tile_subset <- melted_data[Level_ID == level_var, ]
  colnames(tile_subset)[colnames(tile_subset) == 'ID'] <- sample_var
  setkeyv(tile_subset, sample_var)
  setkeyv(metadata, sample_var)
  
  tile_subset <- metadata[tile_subset]
  tile_subset <- tile_subset[!is.na(tile_subset[[group_var]]), ]
  
  if(length(analysis_subset) > 0) {
    tile_subset <- data_subset_long(tile_subset, analysis_subset)
  }
  
  
  sample_order <- unique(tile_subset[order(group_var), sample_var,with=FALSE]) #Error in unique(tile_subset[order(group_var), sample_var]) : error in evaluating the argument 'x' in selecting a method for function 'unique': Error in `[.data.table`(tile_subset, order(group_var), sample_var) : j (the 2nd argument inside [...]) is a single symbol but column name 'sample_var' is not found. Perhaps you intended DT[,..sample_var] or DT[,sample_var,with=FALSE]. This difference to data.frame is deliberate and explained in FAQ 1.1.
  tile_subset <- within(tile_subset, sample_var
                        <- factor(sample_var,
                                  levels=sample_order,
                                  ordered=T))
  
  setkey(tile_subset, Normalized_Count)
  tile_subset <- tile_subset[, sum(Normalized_Count),
                             by=c(group_var, sample_var, 'Name')]
  names(tile_subset)[length(names(tile_subset))] <- 'Normalized_Count'
  
  numselect <- 20
  tile_names <- heatmap_select_top_counts(tile_subset, group_var,
                                          sample_var, numselect)
  name_count <- length(unique(tile_names$Name))
  
  # while(name_count > 50) {
  #   numselect <- numselect - 1
  #   tile_names <- heatmap_select_top_counts(tile_subset, group_var,
  #                                           sample_var, 1)
  #   name_count <- length(unique(tile_names$Name))
  # }
  
  
  tile_subset <- tile_subset[Name %in% tile_names$Name, ]
  
  ## Check if needing to change factor order. Need to convert to data frame first, then back to data table
  tile_subset <- as.data.frame(tile_subset)
  ifelse(factor_order != '', tile_subset[[group_var]] <- factor(tile_subset[[group_var]], levels = eval(factor_order)),1)
  tile_subset <- as.data.table(tile_subset)
  
  # Plot object
  tile <- ggplot(tile_subset, aes_string(x=sample_var, y='Name')) +
    geom_tile(aes(fill=log2(Normalized_Count+1))) +
    facet_wrap(as.formula(paste('~', group_var)), strip.position ='bottom', scales = 'free_x', nrow = 1) +
    theme(strip.text.x=element_text(size=20),
          axis.text.y=element_text(size=20),
          axis.text.x=element_blank(),
          axis.title.x=element_text(size=20),
          axis.title.y=element_blank(),
          legend.position="bottom",
          legend.title=element_text(size=20),
          plot.title=element_text(size=24, hjust=0.5),
          panel.background = element_rect(fill = "white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing=unit(0.1, "lines"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    xlab(paste('Samples by ', group_var, sep='', collapse='')) +
    scale_fill_gradientn(colors=heatmap_col) +
    labs(fill= 'Log2 Normalized Count') +
    ggtitle(paste(data_type, ' ', level_var, ' Normalized Counts by ', group_var, '\n',
                  sep='', collapse=''))
  png(filename=paste(outdir, '/', level_var, '_', group_var, '_',
                     'Heatmap.png', sep='', collapse=''), width=1400, height=700)
  print(tile)
  dev.off()
}

##############################
###                        ###
####   Alpha diversity    ####
###                        ###
##############################

alpha_CSS_diversity <- function(X, minlevel, method='invsimpson') {
  S <- specnumber(X, MARGIN=2)
  alphadiv <- diversity(X, index=method, MARGIN=2)
  return(list(CSS_species_abundance=S,
              CSS_data=X,
              alphadiv=alphadiv))
}


meg_alpha_diversity <- function(data_list,
                                data_names,
                                metadata,
                                sample_var,
                                group_var,
                                analysis_subset,
                                outdir,
                                data_type,
                                factor_order) {
  #Create tables that hold the data
  all_alphadiv <- data.table(ID=character(),
                             Level=character(),
                             Value=numeric())
  names(all_alphadiv)[1] <- sample_var
  
  # CSS_species_abundance
  all_species_CSS <- data.table(ID=character(),
                                Level=character(),
                                Value=numeric())
  names(all_species_CSS)[1] <- sample_var
  
  local_data <- data_list
  
  # Output Diversity counts
  all_sample_data <- alpha_CSS_diversity(MRcounts(local_data[[1]]))
  write.table(all_sample_data$alphadiv, paste(outdir, '/', data_type, '_CSS_InvSimpson_values', group_var, '.csv',
                                              sep='', collapse=''), sep=",")
  
  # I don't currently use the non_zero_sample variable for CSS counts, but it could come in handy for removing those samples or including a report for which samples had 0 counts
  for( l in 1:length(local_data) ) {
    
    if(length(analysis_subset) > 0) {
      local_data[[l]] <- data_subset(local_data[[l]], analysis_subset)
    }
    
    sample_counts <- colSums(MRcounts(local_data[[l]]))
    if(min(sample_counts) == 0) {
      non_zero_sample <- min(sample_counts[sample_counts > 0])
    }
    else {
      non_zero_sample <- 0
    }
    
    # Create main object with alpha rarefaction values.
    # This object includes CSS_species_abundance, CSS_data, and alphadiv
    local_obj <- alpha_CSS_diversity(MRcounts(local_data[[l]]))
    
    # Alpha rarefaction
    temp <- data.table(ID=names(local_obj$alphadiv),
                       Level=rep(data_names[l],
                                 length(local_obj$alphadiv)),
                       Value=as.numeric(local_obj$alphadiv))
    names(temp)[1] <- sample_var
    all_alphadiv <- rbind(all_alphadiv, temp)
    
    # Raw species
    temp <- data.table(ID=names(local_obj$CSS_species_abundance),
                       Level=rep(data_names[l],
                                 length(local_obj$CSS_species_abundance)),
                       Value=as.numeric(local_obj$CSS_species_abundance))
    names(temp)[1] <- sample_var
    all_species_CSS <- rbind(all_species_CSS, temp)
    
  }
  
  #
  ## Alpha rarefaction
  #
  all_alphadiv <- within(all_alphadiv, Level
                         <- factor(Level, levels=data_names,
                                   ordered=T))
  
  setkeyv(all_alphadiv, sample_var)
  setkeyv(metadata, sample_var)
  all_alphadiv <- metadata[all_alphadiv]
  all_alphadiv <- all_alphadiv[!is.na(all_alphadiv[[group_var]]), ]
  
  
  alphadiv_type_sums <- all_alphadiv[Level==data_names[2], median(round(Value, digits=0)), by=group_var]
  alphadiv_value_labels <- as.character(alphadiv_type_sums[[group_var]][order(alphadiv_type_sums$V1, decreasing=T)])
  all_alphadiv[[group_var]] <- factor(all_alphadiv[[group_var]],
                                      levels=alphadiv_value_labels, ordered=T)
  # Check if needing to change factor order.
  ifelse(factor_order != '', all_alphadiv[[group_var]] <- factor(all_alphadiv[[group_var]], levels = eval(factor_order), ordered=T),1)
  
  #print(all_alphadiv)
  #all_alphadiv[[group_var]] <- droplevels( all_alphadiv[[group_var]])
  #all_alphadiv$Level <- droplevels( all_alphadiv$Level)
  # Plot figure
  png(filename=paste(outdir, '/', data_type, '_CSS_alphadiversity_by_', group_var, '.png',
                     sep='', collapse=''),
      width=1024, height=768)
  g_alphadiv <- ggplot(data=all_alphadiv, aes_string(x=group_var,
                                                     y='Value',
                                                     fill=group_var)) +
    geom_boxplot(size=1) +
    facet_wrap(~Level, scales='free_y')
  #g_alphadiv <- g_alphadiv + scale_fill_tableau("Classic 20", direction = -1)
  g_alphadiv <- g_alphadiv +
    ggtitle(paste('Alpha Diversity by ', group_var, ' for CSS data\nInverse Simpson Index',
                  sep='', collapse='')) +
    ylab('Inverse Simpson\'s Index\n') +
    xlab(paste('\n', group_var, sep='', collapse='')) +
    theme(strip.text.x=element_text(size=26),
          axis.text.y=element_text(size=20),
          axis.text.x=element_blank(),
          axis.title.x=element_text(size=26),
          axis.title.y=element_text(size=26, vjust=1),
          #legend.position="right",
          legend.title=element_text(size=24, vjust=1),
          legend.text=element_text(size=20),
          plot.title=element_text(size=30, hjust=0.5),
          panel.background = element_rect(fill = "white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing=unit(0.1, "lines"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))
  print(g_alphadiv)
  dev.off()
  
  
  #
  ## Raw species abundance
  #
  all_species_CSS <- within(all_species_CSS, Level
                            <- factor(Level, levels=data_names,
                                      ordered=T))
  setkeyv(all_species_CSS, sample_var)
  all_species_CSS <- metadata[all_species_CSS]
  all_species_CSS <- all_species_CSS[!is.na(all_species_CSS[[group_var]]), ]
  species_raw_type_sums <- all_species_CSS[Level==data_names[2], median(Value), by=group_var]
  species_raw_value_labels <- as.character(species_raw_type_sums[[group_var]][order(species_raw_type_sums$V1,
                                                                                    decreasing=T)])
  all_species_CSS[[group_var]] <- factor(all_species_CSS[[group_var]],
                                         levels=species_raw_value_labels, ordered=T)
  # Check if needing to change factor order.
  ifelse(factor_order != '', all_species_CSS[[group_var]] <- factor(all_species_CSS[[group_var]], levels = eval(factor_order), ordered=T),1)
  
  # Plot figure
  png(filename=paste(outdir, '/', data_type, '_CSS_richness_by_', group_var, '.png',
                     sep='', collapse=''),
      width=1024, height=768)
  g_sraw <- ggplot(data=all_species_CSS, aes_string(group_var, 'Value', fill=group_var)) +
    geom_boxplot(size=1) +
    facet_wrap(~Level, scales='free_y')
  #g_sraw <- g_sraw + scale_fill_tableau("Classic 20", direction = -1)
  g_sraw <- g_sraw +
    ggtitle(paste('Species Richness by ', group_var, ' for CSS data\nInverse Simpson Index',
                  sep='', collapse='')) +
    ylab('Unique Species\n') +
    xlab(paste('\n', group_var, sep='', collapse='')) +
    theme(strip.text.x=element_text(size=26),
          axis.text.y=element_text(size=20),
          axis.text.x=element_blank(),
          axis.title.x=element_text(size=26),
          axis.title.y=element_text(size=26, vjust=1),
          #legend.position="right",
          legend.title=element_text(size=24, vjust=1),
          legend.text=element_text(size=20),
          plot.title=element_text(size=30, hjust=0.5),
          panel.background = element_rect(fill = "white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing=unit(0.1, "lines"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))
  print(g_sraw)
  dev.off()
  
}






##############################
###                        ###
####   Alpha rarefaction  ####
###                        ###
##############################



# Function that returns species count, rarefied species count, and alpha diversity measures
# for each sample in the m x n matrix, m = features, n = samples
alpha_rarefaction <- function(X, minlevel, method='invsimpson') {
  S <- specnumber(X, MARGIN=2)
  raremax <- min(colSums(X))
  if( raremax < minlevel ) raremax <- minlevel
  Srare <- rarefy(X, raremax, MARGIN=2)
  Xrare <- t(rrarefy(t(X), raremax))
  alphadiv <- diversity(Xrare, index=method, MARGIN=2)
  return(list(raw_species_abundance=S,
              rarefied_species_abundance=Srare,
              rarefied_data=Xrare,
              alphadiv=alphadiv))
}

meg_alpha_rarefaction <- function(data_list,
                                  data_names,
                                  metadata,
                                  sample_var,
                                  group_var,
                                  analysis_subset,
                                  outdir,
                                  data_type,
                                  factor_order) {
  all_alphadiv <- data.table(ID=character(),
                             Level=character(),
                             Value=numeric())
  names(all_alphadiv)[1] <- sample_var
  setkeyv(metadata, sample_var)
  all_species_raw <- data.table(ID=character(),
                                Level=character(),
                                Value=numeric())
  names(all_species_raw)[1] <- sample_var
  
  all_species_rare <- data.table(ID=character(),
                                 Level=character(),
                                 Value=numeric())
  names(all_species_rare)[1] <- sample_var
  
  local_data <- data_list
  
  for( l in 1:length(local_data) ) {
    
    if(length(analysis_subset) > 0) {
      local_data[[l]] <- data_subset(local_data[[l]], analysis_subset)
    }
    
    sample_counts <- colSums(MRcounts(local_data[[l]]))
    if(min(sample_counts) == 0) {
      non_zero_sample <- min(sample_counts[sample_counts > 0])
    }
    else {
      non_zero_sample <- 0
    }
    
    # Create main object with alpha rarefaction values.
    # This object includes raw_species_abundance, rarefied_species_abundance, rarefied_data, and alphadiv
    local_obj <- alpha_rarefaction(MRcounts(local_data[[l]]), minlevel = non_zero_sample)
    
    # Alpha rarefaction
    temp <- data.table(ID=names(local_obj$alphadiv),
                       Level=rep(data_names[l],
                                 length(local_obj$alphadiv)),
                       Value=as.numeric(local_obj$alphadiv))
    names(temp)[1] <- sample_var
    all_alphadiv <- rbind(all_alphadiv, temp)
    
    # Raw species
    temp <- data.table(ID=names(local_obj$raw_species_abundance),
                       Level=rep(data_names[l],
                                 length(local_obj$raw_species_abundance)),
                       Value=as.numeric(local_obj$raw_species_abundance))
    names(temp)[1] <- sample_var
    all_species_raw <- rbind(all_species_raw, temp)
    
    # Rarefied species abundance
    temp <- data.table(ID=names(local_obj$rarefied_species_abundance),
                       Level=rep(data_names[l],
                                 length(local_obj$rarefied_species_abundance)),
                       Value=as.numeric(local_obj$rarefied_species_abundance))
    names(temp)[1] <- sample_var
    all_species_rare <- rbind(all_species_rare, temp)
  }
  
  #
  ## Alpha rarefaction
  #
  all_alphadiv <- within(all_alphadiv, Level
                         <- factor(Level, levels=data_names,
                                   ordered=T))
  setkeyv(all_alphadiv, sample_var)
  all_alphadiv <- metadata[all_alphadiv]
  all_alphadiv <- all_alphadiv[!is.na(all_alphadiv[[group_var]]), ]
  
  alphadiv_type_sums <- all_alphadiv[Level==data_names[2], median(Value), by=group_var]
  alphadiv_value_labels <- as.character(alphadiv_type_sums[[group_var]][order(alphadiv_type_sums$V1, decreasing=T)])
  all_alphadiv[[group_var]] <- factor(all_alphadiv[[group_var]],
                                      levels=alphadiv_value_labels, ordered=T)
  # Check if needing to change factor order.
  ifelse(factor_order != '', all_alphadiv[[group_var]] <- factor(all_alphadiv[[group_var]], levels = eval(factor_order), ordered=T),1)
  
  # Plot figure
  png(filename=paste(outdir, '/', data_type, '_alphadiversity_by_', group_var, '.png',
                     sep='', collapse=''),
      width=1024, height=768)
  g_alphadiv <- ggplot(data=all_alphadiv, aes_string(x=group_var,
                                                     y='Value',
                                                     fill=group_var)) +
    geom_boxplot(size=1) +
    facet_wrap(~Level, scales='free_y')
  #g_alphadiv <- g_alphadiv + scale_fill_tableau("Classic 20", direction = -1)
  g_alphadiv <- g_alphadiv +
    ggtitle(paste('Alpha Diversity by ', group_var, ' for Rarefied data\nInverse Simpson Index',
                  sep='', collapse='')) +
    ylab('Inverse Simpson\'s Index\n') +
    xlab(paste('\n', group_var, sep='', collapse='')) +
    theme(strip.text.x=element_text(size=26),
          axis.text.y=element_text(size=20),
          axis.text.x=element_blank(),
          axis.title.x=element_text(size=26),
          axis.title.y=element_text(size=26, vjust=1),
          #legend.position="right",
          legend.title=element_text(size=24, vjust=1),
          legend.text=element_text(size=20),
          plot.title=element_text(size=30, hjust=0.5),
          panel.background = element_rect(fill = "white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing=unit(0.1, "lines"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))
  print(g_alphadiv)
  dev.off()
  
  
  #
  ## Raw species abundance
  #
  all_species_raw <- within(all_species_raw, Level
                            <- factor(Level, levels=data_names,
                                      ordered=T))
  setkeyv(all_species_raw, sample_var)
  all_species_raw <- metadata[all_species_raw]
  all_species_raw <- all_species_raw[!is.na(all_species_raw[[group_var]]), ]
  species_raw_type_sums <- all_species_raw[Level==data_names[2], median(Value), by=group_var]
  species_raw_value_labels <- as.character(species_raw_type_sums[[group_var]][order(species_raw_type_sums$V1,
                                                                                    decreasing=T)])
  all_species_raw[[group_var]] <- factor(all_species_raw[[group_var]],
                                         levels=species_raw_value_labels, ordered=T)
  # Check if needing to change factor order.
  ifelse(factor_order != '', all_species_raw[[group_var]] <- factor(all_species_raw[[group_var]], levels = eval(factor_order), ordered=T),1)
  
  # Plot figure
  png(filename=paste(outdir, '/', data_type, '_raw_richness_by_', group_var, '.png',
                     sep='', collapse=''),
      width=1024, height=768)
  g_sraw <- ggplot(data=all_species_raw, aes_string(group_var, 'Value', fill=group_var)) +
    geom_boxplot(size=1) +
    facet_wrap(~Level, scales='free_y')
  #g_sraw <- g_sraw + scale_fill_tableau("Classic 20", direction = -1)
  g_sraw <- g_sraw +
    ggtitle(paste('Species Richness by ', group_var, ' for Raw data\nInverse Simpson Index',
                  sep='', collapse='')) +
    ylab('Unique Species\n') +
    xlab(paste('\n', group_var, sep='', collapse='')) +
    theme(strip.text.x=element_text(size=26),
          axis.text.y=element_text(size=20),
          axis.text.x=element_blank(),
          axis.title.x=element_text(size=26),
          axis.title.y=element_text(size=26, vjust=1),
          #legend.position="right",
          legend.title=element_text(size=24, vjust=1),
          legend.text=element_text(size=20),
          plot.title=element_text(size=30, hjust=0.5),
          panel.background = element_rect(fill = "white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing=unit(0.1, "lines"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))
  print(g_sraw)
  dev.off()
  
  
  
  #
  ## Rarefied species abundance
  #
  all_species_rare <- within(all_species_rare, Level
                             <- factor(Level, levels=data_names,
                                       ordered=T))
  setkeyv(all_species_rare, sample_var)
  all_species_rare <- metadata[all_species_rare]
  all_species_rare <- all_species_rare[!is.na(all_species_rare[[group_var]]), ]
  species_rare_type_sums <- all_species_rare[Level==data_names[2], median(Value), by=group_var]
  species_rare_value_labels <- as.character(species_rare_type_sums[[group_var]][order(species_rare_type_sums$V1,
                                                                                      decreasing=T)])
  
  all_species_rare[[group_var]] <- factor(all_species_rare[[group_var]],
                                          levels=species_rare_value_labels, ordered=T)
  ## Check if needing to change factor order.
  ifelse(factor_order != '', all_species_rare[[group_var]] <- factor(all_species_rare[[group_var]], levels = eval(factor_order), ordered=T),1)
  
  png(filename=paste(outdir, '/', data_type, '_rarefied_richness_by_', group_var, '.png',
                     sep='', collapse=''),
      width=1024, height=768)
  g_srare <- ggplot(data=all_species_rare, aes_string(group_var, 'Value', fill=group_var)) +
    geom_boxplot(size=1) +
    facet_wrap(~Level, scales='free_y')
  #g_srare <- g_srare + scale_fill_tableau("Tableau 20", direction = -1)
  g_srare <- g_srare +
    ggtitle(paste('Species Richness by ', group_var, ' for Rarefied data\nInverse Simpson Index',
                  sep='', collapse='')) +
    ylab('Unique Species\n') +
    xlab(paste('\n', group_var, sep='', collapse='')) +
    theme(strip.text.x=element_text(size=26),
          axis.text.y=element_text(size=20),
          axis.text.x=element_blank(),
          axis.title.x=element_text(size=26),
          axis.title.y=element_text(size=26, vjust=1),
          #legend.position="right",
          legend.title=element_text(size=24, vjust=1),
          legend.text=element_text(size=20),
          plot.title=element_text(size=30, hjust=0.5),
          panel.background = element_rect(fill = "white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing=unit(0.1, "lines"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))
  print(g_srare)
  dev.off()
}


##############################
###                        ###
####     Bar plots        ####
###                        ###
##############################
bar_select_top_counts <- function(X, group_var, n) {
  return(X[, tail(.SD, n), by=group_var])
}



meg_relative_barplot <- function(melted_data,
                                 metadata,
                                 sample_var,
                                 group_var,
                                 level_var,
                                 analysis_subset,
                                 outdir,
                                 data_type,
                                 factor_order) {
  setkeyv(melted_data, sample_var)
  setkeyv(metadata, sample_var)
  melted_data <- metadata[melted_data]
  
  if(length(analysis_subset) > 0) {
    bar_subset <- data_subset_long(melted_data, analysis_subset)
  }
  else {
    bar_subset <- melted_data
  }
  
  bar_subset <- data.table(bar_subset[Level_ID == level_var &
                                        !is.na(bar_subset[[group_var]]),
                                      .SD,
                                      .SDcols=c(sample_var,
                                                group_var,
                                                'Name',
                                                'Normalized_Count')])
  setkeyv(bar_subset, sample_var)
  bar_subset <- metadata[bar_subset]
  bar_subset[[sample_var]] <- factor(bar_subset[[sample_var]],
                                     levels=unique(bar_subset[[sample_var]][order(bar_subset[[group_var]])]),
                                     ordered=T)
  
  setkey(bar_subset, Normalized_Count)
  #bar_subset[, sample_number:=(length(unique(bar_subset[[sample_var]]))), by=c(group_var, 'Name')] # commented this one out sample number was the total number of samples
  bar_subset[, sample_number:=(.N), by=c(group_var, 'Name')]
  
  
  bar_subset <- unique(bar_subset[, sum(Normalized_Count) / sample_number,
                                  by=c(group_var, 'Name')])
  ## We should improve how we select the number of taxa to be plotted, for now it's arbitrarily 10 for easy color selection
  numselect <- 21
  bar_names <- bar_select_top_counts(bar_subset, group_var, numselect)
  name_count <- length(unique(bar_names$Name))
  while(name_count > 21) {
    numselect <- numselect - 1
    bar_names <- bar_select_top_counts(bar_subset, group_var, numselect)
    name_count <- length(unique(bar_names$Name))
  }
  
  bar_subset <- bar_subset[Name %in% bar_names$Name, ]
  
  
  names(bar_subset)[length(names(bar_subset))] <- 'Normalized_Count'
  source_sums <- tapply(bar_subset[['Normalized_Count']],
                        bar_subset[[group_var]], sum)
  source_labels <- names(source_sums)[order(source_sums, decreasing=T)]
  bar_subset[[group_var]] <- factor(bar_subset[[group_var]],
                                    levels=source_labels, ordered=T)
  
  # Check if new order of group_var is specified and change if needed
  ifelse(factor_order != '', bar_subset[[group_var]] <- factor(bar_subset[[group_var]], levels = eval(factor_order), ordered=T),1)
  
  # It's not my favorite thing to call the "by" flag with a column, but I couldn't get the group_var variable to work
  bar_subset[, total := sum(Normalized_Count), by = .(bar_subset[[group_var]])][, percentage := Normalized_Count/total , by = .(Name,bar_subset[[group_var]])]
  
  meg_bar <- ggplot(bar_subset, aes_string(x=group_var, y='percentage', fill='Name')) +
    geom_bar(stat='identity') +
    #scale_fill_brewer(palette="Spectral") +
    theme(strip.text.x=element_text(size=18),
          axis.text.y=element_text(size=20),
          axis.text.x=element_text(size=22, vjust=1, hjust=1, angle=33),
          axis.title.x=element_text(size=24),
          axis.title.y=element_text(size=24),
          legend.title=element_text(size=20),
          legend.text=element_text(size=18),
          plot.title=element_text(size=26, hjust=0.25),
          panel.background = element_rect(fill = "white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing=unit(0.1, "lines"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    xlab(group_var) +
    ylab('Relative abundance\n') +
    #scale_fill_tableau("Classic 20", direction = -1) +
    ggtitle(paste(data_type, ' ', level_var, ' Relative Abundance by ', group_var, '\n',
                  sep='', collapse=''))
  png(filename=paste(outdir, '/', data_type, '_', level_var, '_RelativeAbundance_BarPlot_by_', group_var, '.png',
                     sep='', collapse=''), width=1024, height=768)
  print(meg_bar)
  dev.off()
}






meg_barplot <- function(melted_data,
                        metadata,
                        sample_var,
                        group_var,
                        level_var,
                        analysis_subset,
                        outdir,
                        data_type,
                        factor_order) {
  setkeyv(melted_data, sample_var)
  setkeyv(metadata, sample_var)
  melted_data <- metadata[melted_data]
  
  if(length(analysis_subset) > 0) {
    bar_subset <- data_subset_long(melted_data, analysis_subset)
  }
  else {
    bar_subset <- melted_data
  }
  
  bar_subset <- data.table(bar_subset[Level_ID == level_var &
                                        !is.na(bar_subset[[group_var]]),
                                      .SD,
                                      .SDcols=c(sample_var,
                                                group_var,
                                                'Name',
                                                'Normalized_Count')])
  setkeyv(bar_subset, sample_var)
  bar_subset <- metadata[bar_subset]
  bar_subset[[sample_var]] <- factor(bar_subset[[sample_var]],
                                     levels=unique(bar_subset[[sample_var]][order(bar_subset[[group_var]])]),
                                     ordered=T)
  
  setkey(bar_subset, Normalized_Count)
  #bar_subset[, sample_number:=(length(unique(bar_subset[[sample_var]]))), by=c(group_var, 'Name')] # commented this one out, sample number was the total number of samples
  bar_subset[, sample_number:=(.N), by=c(group_var, 'Name')]
  
  
  bar_subset <- unique(bar_subset[, sum(Normalized_Count) / sample_number,
                                  by=c(group_var, 'Name')])
  
  numselect <- 11
  bar_names <- bar_select_top_counts(bar_subset, group_var, numselect)
  name_count <- length(unique(bar_names$Name))
  while(name_count > 11) {
    numselect <- numselect - 1
    bar_names <- bar_select_top_counts(bar_subset, group_var, numselect)
    name_count <- length(unique(bar_names$Name))
  }
  
  bar_subset <- bar_subset[Name %in% bar_names$Name, ]
  
  
  names(bar_subset)[length(names(bar_subset))] <- 'Normalized_Count'
  source_sums <- tapply(bar_subset[['Normalized_Count']],
                        bar_subset[[group_var]], sum)
  source_labels <- names(source_sums)[order(source_sums, decreasing=T)]
  bar_subset[[group_var]] <- factor(bar_subset[[group_var]],
                                    levels=source_labels, ordered=T)
  
  # Check if new order of group_var is specified and change if needed
  ifelse(factor_order != '', bar_subset[[group_var]] <- factor(bar_subset[[group_var]], levels = eval(factor_order), ordered=T),1)
  
  bar_subset$Name <-  droplevels(bar_subset$Name)
  
  meg_bar <- ggplot(bar_subset, aes_string(x=group_var, y='Normalized_Count', fill='Name')) +
    geom_bar(stat='identity') +
    scale_fill_brewer(palette="Spectral") +
    theme(strip.text.x=element_text(size=18),
          axis.text.y=element_text(size=20),
          axis.text.x=element_text(size=22, vjust=1, hjust=1, angle=33),
          axis.title.x=element_text(size=24),
          axis.title.y=element_text(size=24),
          legend.title=element_text(size=20),
          legend.text=element_text(size=18),
          plot.title=element_text(size=26, hjust=0.25),
          panel.background = element_rect(fill = "white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing=unit(0.1, "lines"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    xlab(group_var) +
    ylab('Mean of Normalized Count\n') +
    ggtitle(paste('Mean ', data_type, ' ', level_var, ' Normalized Count by ', group_var, '\n',
                  sep='', collapse=''))
  png(filename=paste(outdir, '/', data_type, '_', level_var, '_Mean_BarPlot_by_', group_var, '.png',
                     sep='', collapse=''), width=1024, height=768)
  print(meg_bar)
  dev.off()
}





meg_median_barplot <- function(melted_data,
                               metadata,
                               sample_var,
                               group_var,
                               level_var,
                               analysis_subset,
                               outdir,
                               data_type,
                               factor_order) {
  setkeyv(melted_data, sample_var)
  setkeyv(metadata, sample_var)
  melted_data <- metadata[melted_data]
  
  if(length(analysis_subset) > 0) {
    bar_subset <- data_subset_long(melted_data, analysis_subset)
  }
  else {
    bar_subset <- melted_data
  }
  
  bar_subset <- data.table(bar_subset[Level_ID == level_var &
                                        !is.na(bar_subset[[group_var]]),
                                      .SD,
                                      .SDcols=c(sample_var,
                                                group_var,
                                                'Name',
                                                'Normalized_Count')])
  setkeyv(bar_subset, sample_var)
  bar_subset <- metadata[bar_subset]
  bar_subset[[sample_var]] <- factor(bar_subset[[sample_var]],
                                     levels=unique(bar_subset[[sample_var]][order(bar_subset[[group_var]])]),
                                     ordered=T)
  
  setkey(bar_subset, Normalized_Count)
  #bar_subset[, sample_number:=(length(unique(bar_subset[[sample_var]]))), by=c(group_var, 'Name')] # commented this one out, sample number was the total number of samples
  bar_subset[, sample_number:=(.N), by=c(group_var, 'Name')]
  
  
  bar_subset <- unique(bar_subset[, median(round(Normalized_Count)),
                                  by=c(group_var, 'Name')])
  
  numselect <- 11
  bar_names <- bar_select_top_counts(bar_subset, group_var, numselect)
  name_count <- length(unique(bar_names$Name))
  while(name_count > 11) {
    numselect <- numselect - 1
    bar_names <- bar_select_top_counts(bar_subset, group_var, numselect)
    name_count <- length(unique(bar_names$Name))
  }
  
  bar_subset <- bar_subset[Name %in% bar_names$Name, ]
  
  
  names(bar_subset)[length(names(bar_subset))] <- 'Normalized_Count'
  source_sums <- tapply(bar_subset[['Normalized_Count']],
                        bar_subset[[group_var]], sum)
  source_labels <- names(source_sums)[order(source_sums, decreasing=T)]
  bar_subset[[group_var]] <- factor(bar_subset[[group_var]],
                                    levels=source_labels, ordered=T)
  
  # Check if new order of group_var is specified and change if needed
  ifelse(factor_order != '', bar_subset[[group_var]] <- factor(bar_subset[[group_var]], levels = eval(factor_order), ordered=T),1)
  
  bar_subset$Name <-  droplevels(bar_subset$Name)
  
  meg_bar <- ggplot(bar_subset, aes_string(x=group_var, y='Normalized_Count', fill='Name')) +
    geom_bar(stat='identity') +
    scale_fill_brewer(palette="Spectral") +
    theme(strip.text.x=element_text(size=18),
          axis.text.y=element_text(size=20),
          axis.text.x=element_text(size=22, vjust=1, hjust=1, angle=33),
          axis.title.x=element_text(size=24),
          axis.title.y=element_text(size=24),
          legend.title=element_text(size=20),
          legend.text=element_text(size=18),
          plot.title=element_text(size=26, hjust=0.25),
          panel.background = element_rect(fill = "white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing=unit(0.1, "lines"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    xlab(group_var) +
    ylab('Normalized Count\n') +
    ggtitle(paste('Median ', data_type, ' ', level_var, ' Normalized Count by ', group_var, '\n',
                  sep='', collapse=''))
  png(filename=paste(outdir, '/', data_type, '_', level_var, '_Median_BarPlot_by_', group_var, '.png',
                     sep='', collapse=''), width=1024, height=768)
  print(meg_bar)
  dev.off()
}

##############################
###                        ###
####     ZIG model        ####
###                        ###
##############################

meg_fitZig <- function(data_list,
                       data_names,
                       metadata,
                       zero_mod,
                       data_mod,
                       filter_min_threshold,
                       contrast_list,
                       random_effect_var,
                       outdir,
                       analysis_name,
                       analysis_subset,
                       data_type,
                       pval=0.1,
                       top_hits=200) {
  settings <- zigControl(maxit=20000, verbose=F)
  
  local_obj <- data_list
  res <- list()
  for( l in 1:length(local_obj) ) {
    filter_threshold <- quantile(rowSums(MRcounts(local_obj[[l]])), 0.15)
    if( filter_threshold > filter_min_threshold ) filter_threshold <- filter_min_threshold
    local_obj[[l]] <- local_obj[[l]][which(rowSums(MRcounts(local_obj[[l]])) >= filter_threshold ), ]
    if(length(analysis_subset) > 0) {
      local_obj[[l]] <- data_subset(local_obj[[l]], analysis_subset)
      amr_sparseFeatures = which(rowSums(MRcounts(local_obj[[l]]) > 0) < 1) ## Check to see if any features in the subset == 0
      if (length(amr_sparseFeatures) != 0) {
        local_obj[[l]] = local_obj[[l]][-amr_sparseFeatures,]
      }
    }
    
    
    col_selection <- as.integer(which(colSums(MRcounts(local_obj[[l]]) > 0) > 1))
    local_obj[[l]] <- local_obj[[l]][, col_selection]
    
    mod_select <- model.matrix(eval(parse(text=data_mod)), data=pData(local_obj[[l]]))
    zero_mod_select <- zero_mod[col_selection, ]
    
    cumNorm(local_obj[[l]])  # This is a placeholder for metagenomeSeq; we don't actually use these values
    
    tryCatch(
      {
        if( is.na(random_effect_var) ) {
          res[[l]] <- fitZig(obj=local_obj[[l]],
                             mod=mod_select,
                             zeroMod=zero_mod_select,
                             control=settings,
                             useCSSoffset=F)
        }
        else {
          res[[l]] <- fitZig(obj=local_obj[[l]],
                             mod=mod_select,
                             zeroMod=zero_mod_select,
                             control=settings,
                             useCSSoffset=F,
                             useMixedModel=T,
                             block=pData(local_obj[[l]])[, random_effect_var])
        }
      },
      error=function(e) {
        print(paste('Model failed to converge for ', data_type, ' ', data_names[l], ' ', analysis_name,
                    sep='', collapse=''))
      },
      finally={
        if( length(res) != l ) {
          next
        }
      }
    )
    
    local_contrasts <- contrast_list
    local_contrasts[[length(local_contrasts)+1]] <- res[[l]]$fit$design
    names(local_contrasts)[length(local_contrasts)] <- 'levels'
    
    contrast_matrix <- do.call(makeContrasts, local_contrasts)
    colnames(contrast_matrix) <- make.names(contrast_list)
    
    contrast_fit <- contrasts.fit(res[[l]]$fit, contrast_matrix)
    contrast_fit <- eBayes(contrast_fit)
    
    stats_results <- data.table(
      Node.Name=character(),
      Contrast=character(),
      logFC=numeric(),
      CI.L=numeric(),
      CI.R=numeric(),
      AveExpr=numeric(),
      t=numeric(),
      P.Value=numeric(),
      adj.P.Val=numeric(),
      B=numeric()
    )
    
    for( c in 1:ncol(contrast_fit$contrasts) ) {
      tophits <- topTable(contrast_fit, p.value=pval, confint=T,
                          number=top_hits, sort.by='AveExpr', coef=c)
      
      if( nrow(tophits) > 0) {
        temp_res <- data.table(
          Node.Name=rownames(tophits),
          Contrast=rep(colnames(contrast_fit$contrasts)[c], nrow(tophits))
        )
        temp_res <- cbind(temp_res, tophits)
        stats_results <- rbind(stats_results, temp_res)
      }
      else {
        print(paste('No significant results for', data_type,
                    data_names[l], analysis_name,
                    colnames(contrast_fit$contrasts)[c],
                    sep=' ', collapse=''))
      }
    }
    
    if( nrow(stats_results) > 0 ) {
      write.csv(stats_results,
                file=paste(outdir, '/', analysis_name, '_', data_type, '_',
                           data_names[l], '_',
                           contrast_list[1], '_Model_Contrasts.csv',
                           sep='', collapse=''),
                quote=F, row.names=F)
    }
  }
}



##
### facet_nested.R function from https://github.com/teunbrand/ggh4x/edit/master/R/facet_nested.R
##

# Main function -----------------------------------------------------------

#' @title Layout panels in a grid with nested strips
#'
#' @description \code{facet_nested()} forms a matrix of panels defined by row
#'   and column faceting variables and nests grouped facets.
#'
#' @inheritParams ggplot2::facet_grid
#' @param nest_line a \code{logical} vector of length 1, indicating whether to
#'   draw a nesting line to indicate the nesting of variables. Control the look
#'   of the nesting line by setting the \code{ggh4x.facet.nestline} theme
#'   element.
#' @param resect  a \code{unit} vector of length 1, indicating how much the
#'   nesting line should be shortened.
#' @param bleed a \code{logical} vector of length 1, indicating whether merging
#'   of lower-level variables is allowed when the higher-level variables are
#'   separate. See details.
#'
#' @details Unlike \code{facet_grid()}, this function only automatically expands
#'   missing variables when they have no variables in that direction, to allow
#'   for unnested variables. It still requires at least one layer to have all
#'   faceting variables.
#'
#'   Hierarchies are inferred from the order of variables supplied to
#'   \code{rows} or \code{cols}. The first variable is interpreted to be the
#'   outermost variable, while the last variable is interpreted to be the
#'   innermost variable. They display order is always such that the outermost
#'   variable is placed the furthest away from the panels. Strips are
#'   automatically grouped when they span a nested variable.
#'
#'   The \code{bleed} argument controls whether lower-level variables are allowed
#'   to be merged when higher-level are different, i.e. they can bleed over
#'   hierarchies. Suppose the \code{facet_grid()} behaviour would be the
#'   following:
#'
#'   \code{[_1_][_2_][_2_]} \cr \code{[_3_][_3_][_4_]}
#'
#'   In such case, the default \code{bleed = FALSE} argument would result in the
#'   following:
#'
#'   \code{[_1_][___2____]} \cr \code{[_3_][_3_][_4_]}
#'
#'   Whereas \code{bleed = TRUE} would allow the following:
#'
#'   \code{[_1_][___2____]} \cr \code{[___3____][_4_]}
#'
#' @export
#'
#' @return A \emph{FacetNested} ggproto object.
#' @family facetting functions
#'
#' @seealso See \code{\link[ggplot2]{facet_grid}} for descriptions of the
#'   original arguments. See \code{\link[grid]{unit}} for the construction of a
#'   \code{unit} vector.
#'
#' @examples
#' df <- iris
#' df$nester <- ifelse(df$Species == "setosa",
#'                     "Short Leaves",
#'                     "Long Leaves")
#'
#' ggplot(df, aes(Sepal.Length, Petal.Length)) +
#'   geom_point() +
#'   facet_nested(~ nester + Species)
#'
#' # Controlling the nest line
#' ggplot(df, aes(Sepal.Length, Petal.Length)) +
#'   geom_point() +
#'   facet_nested(~ nester + Species, nest_line = TRUE) +
#'   theme(ggh4x.facet.nestline = element_line(linetype = 3))
facet_nested <- function(
  rows = NULL, cols = NULL, scales = "fixed", space = "fixed",
  shrink = TRUE, labeller = "label_value", as.table = TRUE,
  switch = NULL, drop = TRUE, margins = FALSE, facets = NULL,
  nest_line = FALSE, resect = unit(0, "mm"), bleed = FALSE
) {
  if (!is.null(facets)) {
    rows <- facets
  }
  if (is.logical(cols)) {
    margins <- cols
    cols <- NULL
  }
  scales <- match.arg(scales, c("fixed", "free_x", "free_y", "free"))
  free <- list(x = any(scales %in% c("free_x", "free")),
               y = any(scales %in% c("free_y", "free")))
  
  space <- match.arg(space, c("fixed", "free_x", "free_y", "free"))
  space_free <- list(x = any(space %in% c("free_x", "free")),
                     y = any(space %in% c("free_y", "free")))
  
  if (!is.null(switch) && !switch %in% c("both", "x", "y")) {
    stop("switch must be either 'both', 'x', or 'y'", call. = FALSE)
  }
  
  facets_list <- .int$grid_as_facets_list(rows, cols)
  n <- length(facets_list)
  if (n > 2L) {
    stop("A grid facet specification can't have more than two dimensions",
         .call = FALSE)
  }
  if (n == 1L) {
    rows <- quos()
    cols <- facets_list[[1]]
  } else {
    rows <- facets_list[[1]]
    cols <- facets_list[[2]]
  }
  labeller <- .int$check_labeller(labeller)
  ggproto(NULL, FacetNested, shrink = shrink,
          params = list(
            rows = rows,
            cols = cols,
            margins = margins,
            free = free,
            space_free = space_free,
            labeller = labeller,
            as.table = as.table,
            switch = switch,
            drop = drop,
            nest_line = nest_line,
            resect = resect,
            bleed = bleed
          ))
}

# ggproto -----------------------------------------------------------------

#' @usage NULL
#' @format NULL
#' @export
#' @rdname ggh4x_extensions
FacetNested <- ggproto(
  "FacetNested", FacetGrid,
  map_data = function(data, layout, params) {
    # Handle empty data
    if (.int$empty(data)) {
      return(cbind(data, PANEL = integer(0)))
    }
    # Setup variables
    rows <- params$rows
    cols <- params$cols
    
    vars <- c(names(rows), names(cols))
    if (length(vars) == 0) {
      data$PANEL <- layout$PANEL
      return(data)
    }
    
    margin_vars <- list(intersect(names(rows), names(data)),
                        intersect(names(cols), names(data)))
    
    # Add variables
    data <- .int$reshape_add_margins(data, margin_vars, params$margins)
    facet_vals <- .int$eval_facets(c(rows, cols), data, params$.possible_columns)
    
    # Only set as missing if it has no variable in that direction
    missing_facets <- character(0)
    if (!any(names(rows) %in% names(facet_vals))){
      missing_facets <- c(missing_facets,
                          setdiff(names(rows), names(facet_vals)))
    }
    if (!any(names(cols) %in% names(facet_vals))){
      missing_facets <- c(missing_facets,
                          setdiff(names(cols), names(facet_vals)))
    }
    
    # Fill in missing values
    if (length(missing_facets) > 0) {
      to_add <- unique(layout[missing_facets])
      data_rep <- rep.int(1:nrow(data), nrow(to_add))
      facet_rep <- rep(1:nrow(to_add), each = nrow(data))
      data <- data[data_rep, , drop = FALSE]
      rownames(data) <- NULL
      facet_vals <- cbind(facet_vals[data_rep, , drop = FALSE],
                          to_add[facet_rep, , drop = FALSE])
      rownames(facet_vals) <- NULL
    }
    
    # Match columns to facets
    if (nrow(facet_vals) == 0) {
      data$PANEL <- -1
    } else {
      facet_vals[] <- lapply(facet_vals[], as.factor)
      facet_vals[] <- lapply(facet_vals[], addNA, ifany = TRUE)
      keys <- plyr::join.keys(facet_vals, layout,
                              by = vars[vars %in% names(facet_vals)])
      data$PANEL <- layout$PANEL[match(keys$x, keys$y)]
    }
    data
  },
  compute_layout = function(data, params) {
    rows <- params$rows
    cols <- params$cols
    dups <- intersect(names(rows), names(cols))
    
    if (length(dups) > 0) {
      stop("Facetting variables can only appear in row or cols, not both.\n",
           "Problems: ", paste0(dups, collapse = "'"), call. = FALSE)
    }
    
    base_rows <- combine_nested_vars(data, params$plot_env,
                                     rows, drop = params$drop)
    if (!params$as.table) {
      rev_order <- function(x) factor(x, levels = rev(.int$ulevels(x)))
    }
    base_cols <- combine_nested_vars(data, params$plot_env, cols,
                                     drop = params$drop)
    base <- .int$df.grid(base_rows, base_cols)
    
    if (nrow(base) == 0) {
      return(.int$new_data_frame(list(PANEL = 1L, ROW = 1L, COL = 1L,
                                      SCALE_X = 1L, SCALE_Y = 1L)))
    }
    
    base <- .int$reshape_add_margins(
      base, list(names(rows), names(cols)), params$margins
    )
    base <- unique(base)
    
    panel <- .int$id(base, drop = TRUE)
    panel <- factor(panel, levels = seq_len(attr(panel, "n")))
    
    rows <- if (!length(names(rows))) {
      rep(1L, length(panel))
    } else {
      .int$id(base[names(rows)], drop = TRUE)
    }
    cols <- if (!length(names(cols))) {
      rep(1L, length(panel))
    } else {
      .int$id(base[names(cols)], drop = TRUE)
    }
    
    panels <- .int$new_data_frame(
      c(list(PANEL = panel, ROW = rows, COL = cols), base)
    )
    panels <- panels[order(panels$PANEL), , drop = FALSE]
    rownames(panels) <- NULL
    panels$SCALE_X <- if (params$free$x) {
      panels$COL
    } else {
      1L
    }
    panels$SCALE_Y <- if (params$free$y) {
      panels$ROW
    } else {
      1L
    }
    panels
  },
  draw_panels = function(panels, layout, x_scales, y_scales, ranges, coord,
                         data, theme, params) {
    panel_table <- FacetGrid$draw_panels(panels, layout, x_scales, y_scales,
                                         ranges, coord, data, theme, params)
    
    # Setup strips
    col_vars  <- unique(layout[names(params$cols)])
    row_vars  <- unique(layout[names(params$rows)])
    attr(col_vars, "type")  <- "cols"
    attr(col_vars, "facet") <- "grid"
    attr(row_vars, "type")  <- "rows"
    attr(row_vars, "facet") <- "grid"
    
    # Build strips
    switch_x <- !is.null(params$switch) && params$switch %in% c("both", "x")
    switch_y <- !is.null(params$switch) && params$switch %in% c("both", "y")
    
    # Merging strips
    merge_cols <- apply(col_vars, 2, function(x) any(rle(x)$lengths > 1))
    merge_rows <- apply(row_vars, 2, function(x) any(rle(x)$lengths > 1))
    
    if (any(merge_cols)) {
      if (switch_x) {
        panel_table <- merge_strips(panel_table,
                                    col_vars, switch_x, params, theme, "b")
      } else {
        panel_table <- merge_strips(panel_table,
                                    col_vars, switch_x, params, theme, "t")
      }
    }
    
    if (any(merge_rows)) {
      if (switch_y) {
        panel_table <- merge_strips(panel_table,
                                    row_vars, switch_y, params, theme, "l")
      } else {
        panel_table <- merge_strips(panel_table,
                                    row_vars, switch_y, params, theme, "r")
      }
    }
    panel_table
  }
)

# Helper functions -----------------------------------------------

combine_nested_vars <- function(
  data, env = emptyenv(), vars = NULL, drop = TRUE
) {
  if (length(vars) == 0) {
    return(.int$new_data_frame())
  }
  
  possible_columns <- unique(unlist(lapply(data, names)))
  
  values <- .int$compact(lapply(data, .int$eval_facets, facets = vars,
                                possible_columns = possible_columns))
  has_all <- unlist(lapply(values, length)) == length(vars)
  if (!any(has_all)) {
    missing <- lapply(values, function(x) setdiff(names(vars), names(x)))
    missing_txt <- vapply(missing, .int$var_list, character(1))
    name <- c("Plot", paste0("Layer ", seq_len(length(data) - 1)))
    stop("At least one layer must contain all faceting variables: ",
         .int$var_list(names(vars)), ".\n", paste0("* ", name, " is missing ",
                                                   missing_txt, collapse = "\n"),
         call. = FALSE)
  }
  base <- unique(.int$rbind_dfs(values[has_all]))
  if (!drop) {
    base <- .int$unique_combs(base)
  }
  for (value in values[!has_all]) {
    if (.int$empty(value))
      next
    old <- base[setdiff(names(base), names(value))]
    new <- unique(value[intersect(names(base), names(value))])
    if (drop) {
      new <- .int$unique_combs(new)
    }
    old[setdiff(names(base), names(value))] <- rep("", nrow(old))
    base <- rbind(base, .int$df.grid(old, new))
  }
  if (.int$empty(base)) {
    stop("Facetting variables must have at least one value",
         call. = FALSE)
  }
  base
}

# New merge strips --------------------------------------------------------

merge_strips <- function(
  panel_table, vars, switch, params, theme, where = "t"
) {
  orient <- if (where %in% c("t", "b")) "x" else "y"
  nlevels <- ncol(vars)
  
  these_strips <- grep(paste0("strip-", where), panel_table$layout$name)
  strp_rows <- range(panel_table$layout$t[these_strips])
  strp_cols <- range(panel_table$layout$l[these_strips])
  strp_rows <- seq(strp_rows[1], strp_rows[2])
  strp_cols <- seq(strp_cols[1], strp_cols[2])
  strp <- panel_table[strp_rows, strp_cols]
  
  # Make empty template
  template <- strp
  template$grobs <- list()
  template$layout <- template$layout[0,]
  
  # Inflate strips
  for (i in seq_along(strp$grobs)) {
    sub <- strp$grobs[[i]]
    if (where == "b") {
      sub$layout$t <- rev(sub$layout$t)
      sub$layout$b <- rev(sub$layout$b)
    }
    n <- length(sub$grobs)
    lay <- strp$layout[i,]
    lay <- lay[rep(1, n),]
    rownames(lay) <- NULL
    sub <- lapply(seq_len(n), function(j) {
      x <- sub
      x$grobs <- x$grobs[j]
      x$layout <- x$layout[j,]
      x
    })
    template <- gtable_add_grob(
      template,
      sub, t = lay$t, l = lay$l, b = lay$b, r = lay$r,
      z = lay$z, clip = lay$clip, name = paste0(lay$name, "-", seq_len(n))
    )
  }
  
  if (!params$bleed) {
    vars[] <- lapply(seq_len(ncol(vars)), function(i) {
      do.call(paste0, vars[, seq(i), drop = FALSE])
    })
  }
  merge <- apply(vars, 2, function(x) any(rle(x)$lengths > 1))
  
  if (where == "r") {
    vars <- rev(vars)
    merge <- rev(merge)
  }
  
  # Abstract away strips
  strip_ids <- strsplit(template$layout$name, "-", fixed = TRUE)
  strip_ids <- do.call(rbind, strip_ids)
  strip_ids <- strip_ids[,3:ncol(strip_ids)]
  mode(strip_ids) <- "integer"
  
  template$layout$delete <- rep(FALSE, nrow(strip_ids))
  template$layout$aquire <- seq_along(template$grobs)
  
  for (i in seq_len(nlevels)) {
    if (!merge[i]) {
      next()
    }
    ii <- strip_ids[, 2] == i
    
    # Figure out what to merge
    j <- as.numeric(as.factor(vars[, i]))
    
    ends <- cumsum(rle(j)$lengths)
    starts <- c(1, which(diff(j) != 0) + 1)
    
    # Figure out what strip to remove
    seqs <- unlist(Map(seq, from = starts, to = ends))
    delete_this <- seqs[!(seqs %in% starts)]
    delete_this <- which(strip_ids[, 1] %in% delete_this & ii)
    template$layout$delete[delete_this] <- TRUE
    
    # Figure out what cells to expand
    expand <- seqs[seqs %in% starts]
    expand <- which(strip_ids[, 1] %in% expand & ii)
    expand_where <- seqs[seqs %in% ends]
    expand_where <- which(strip_ids[, 1] %in% expand_where & ii)
    template$layout$aquire[expand] <- template$layout$aquire[expand_where]
  }
  
  # Do expansion
  if (orient == "x") {
    template$layout$r <- template$layout$r[template$layout$aquire]
  } else {
    template$layout$b <- template$layout$b[template$layout$aquire]
  }
  
  # Do deletion
  template$grobs  <- template$grobs[!template$layout$delete]
  strip_ids <- strip_ids[!template$layout$delete,]
  template$layout <- template$layout[!template$layout$delete,]
  
  # Add nesting indicator
  if (params$nest_line) {
    active <- unit(c(0, 1), "npc") + c(1, -1) * params$resect
    passive <- if (switch) c(1, 1) else c(0, 0)
    nindi <- element_render(
      theme, "ggh4x.facet.nestline",
      x = switch(orient, x = active,  y = passive),
      y = switch(orient, x = passive, y = active)
    )
    i <- which(with(template$layout, t != b | l != r))
    offset <- switch(
      orient,
      x = vapply(template$grobs, function(grob){grob$layout$t}, numeric(1)),
      y = vapply(template$grobs, function(grob){grob$layout$l}, numeric(1))
    )
    offset <- if (where %in% c("r", "b")) offset else nlevels - offset
    template$grobs[i] <- lapply(template$grobs[i], function(grb) {
      grb <- with(grb$layout, gtable_add_grob(
        grb, nindi, t = t, l = l, r = r, b = b,
        z = z,
        name = "nester",
        clip = "off"
      ))
    })
    template$layout$z <- template$layout$z + offset
  }
  
  # Delete old strips
  panel_table <- gtable_filter(panel_table, paste0("strip-", where),
                               fixed = TRUE, trim = FALSE, invert = TRUE)
  
  # Place back new strips
  panel_table <- with(template$layout, gtable_add_grob(
    panel_table,
    template$grobs,
    t = t - 1 + strp_rows[1],
    l = l - 1 + strp_cols[1],
    b = b - 1 + strp_rows[1],
    r = r - 1 + strp_cols[1],
    z = z, clip = clip, name = name
  ))
  panel_table
}
