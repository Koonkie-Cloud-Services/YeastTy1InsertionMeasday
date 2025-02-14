
extract_metadata_from_library_names <- function(df) {
  df %>% 
    mutate(
      subunit = case_when(
        grepl("34", library) ~ "RPC34",
        grepl("53", library) ~ "RPC53",
        grepl("160", library) ~ "RPC160",
        .default = "Other"),
      sugar = case_when(
        grepl(x=library, pattern = "_Dex$", perl = TRUE, ignore.case = TRUE) ~ "Dextrose",
        grepl(x=library, pattern = "_Gal$", perl = TRUE, ignore.case = TRUE) ~ "Galactose",
        .default = NA
      ),
      strain.type = case_when(
        grepl(x=library, pattern = "WT_", perl = TRUE, ignore.case = TRUE) ~ "Wild-type",
        .default = "Mutant"
      ),
      group = gsub(x = library, pattern = "_(Gal|Dex)$", replacement = ''),
      group = gsub(x = group, pattern = "(?i)r[op]c[0-9]{2,3}", replacement = ''),
      group = stringr::str_remove(group, pattern = "deletion_set_background"),
      group = gsub(x = group, pattern = "^(\\.|_)", replacement = ''),
      group = gsub(x = group, pattern = "\\.", replacement = '_'),
    ) %>% 
    return()
}

.load_genome_coverage_table <- function(file_path, col_name_lookup) {
  read.table(file = file_path,
             header = FALSE,
             row.names = NULL,
             sep = "\t",
             comment.char = "",
             quote = "") %>%
    dplyr::rename(tidyselect::any_of(col_name_lookup)) %>% 
    dplyr::filter(chr != "chrmt") %>% 
    return()
}

load_genome_coverage_data_frame <- function(file_path) {
  hm.cov_file_column_names <- c(
    "library"="V1",
    "chr"="V2",
    "start"="V3",
    "end"="V4",
    "feature.name"="V5",
    "strand"="V6",
    "coverage"="V7"
  )
  
  .load_genome_coverage_table(
    file_path = file_path,
    col_name_lookup = hm.cov_file_column_names
  ) %>%
    dplyr::mutate(feature.length = end - start) %>% 
    .fix_gene_names() %>%
    filter_blackout_regions(coord_df = blackout_df) %>%
    # dplyr::filter(feature.length > 1) %>%  # No longer necessary to remove anti-features
    calculate_tpm(read_counts %>% 
                    dplyr::filter(norm_type == p.read_norm_type)) %>% 
    dplyr::group_by(library, feature.name) %>% 
    dplyr::reframe(
      read_count = sum(coverage),
      sum_cov = sum(coverage.pbp),
      mean_cov = mean(coverage.pbp),
      mean_cov.norm = mean(coverage.pbp.norm),
      sd_cov.norm = sd(coverage.pbp.norm),
      mean_tpm = mean(tpm)
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange(library, feature.name) %>% 
    return()
}

calculate_tpm <- function(cov_df,
                          lib_size_df,
                          drop_unaligned_regions = TRUE) {
  cov_df <- .fix_library_names(cov_df)
  lib_size_df <- .fix_library_names(lib_size_df)
  
  if (length(dplyr::setdiff(cov_df$library, lib_size_df$library)) > 0) {
    stop(
      paste("Library names in cov_df and lib_size_df do not match:\n", 
            paste(dplyr::setdiff(cov_df$library, lib_size_df$library),
                  collapse = ','))
    )
  }
  if (drop_unaligned_regions == TRUE) {
    cov_df <- dplyr::filter(cov_df, coverage > 0)
  }
  
  dplyr::left_join(x = cov_df,
                   y = lib_size_df,
                   by = "library") %>%
    dplyr::ungroup() %>% 
    dplyr::mutate(
      # coverage = coverage - 1,
      rpk = coverage / (feature.length / 1E3)) %>% 
    dplyr::group_by(library) %>% 
    dplyr::mutate(scaler = sum(rpk) / 1E6) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(
      tpm = rpk / scaler,
      coverage.pbp = coverage / feature.length,
      coverage.pbp.norm = coverage.pbp / ( seq_library_size / 1E6)
    ) %>%
    dplyr::select(!tidyselect::all_of(
      c("rpk", "scaler", "norm_type", "seq_library_size")
    )) %>% 
  return()
}

.fix_gene_names <- function(df) {
  df %>% 
    dplyr::mutate(
      gene.name = ifelse(grepl(pattern = "^.*(gene))",
                               x = feature.name),
                         stringr::str_extract(string = feature.name,
                                              pattern = "(.*)\\(gene\\)",
                                              group = 1),
                         NA),
      feature.name = gsub(pattern = "^.*(gene))",
                          replacement = "gene",
                          x = feature.name),
      feature.name = gsub("MRibosomal",
                          "Ribosomal",
                          feature.name),
      feature.name = gsub("_", " ", feature.name),
      feature.name = case_when(
        feature.name == "Ribosomal Protein" ~ "RPL",
        feature.name == "intergenic" ~ "anti-feature",
        feature.name == "rep origin" ~ "origin of replication",
        .default = feature.name
      )
    ) %>%
    return()
}

.fix_library_names <- function(df) {
  df %>% dplyr::mutate(
    library = gsub("&|\\.|-",
                   "_",
                   library)
  ) %>% 
    return()
}

.subtract_antifeature_from_feature_freqs <- function(df) {
  # TODO: implement if necessary
  # all_coverage_info_agg_dcast_df[-1, ] = all_coverage_info_agg_dcast_df[-1 , ] - all_coverage_info_agg_dcast_df[rep(1, nrow(all_coverage_info_agg_dcast_df) - 1), ]
  # all_coverage_info_agg_dcast_df[1,] = 0
  df %>% 
    return()
}

filter_blackout_regions <- function(df, coord_df) {
  for (i in nrow(coord_df)) {
    r <- coord_df[i,]
    df <- df %>% 
      filter(chr != r$chr | end < r$start | start > r$end)
  }
  return(df)
}

colour_metadata <- function(df,
                            palette_vals = NULL,
                            brewer_palette_name = "Dark2",
                            palette_size = NULL) {
  
  # Select the columns with discrete values
  df <- df %>%
    as.data.frame() %>%
    dplyr::select(!dplyr::where(is.numeric))
  if (ncol(df) == 0) {
    return(NA)
  }
  
  if ( length(brewer_palette_name) > 1 ) {
    if ( length(brewer_palette_name) != ncol(df) ) {
      stop(paste(
        "Number of palette names must equal number of categorical columns in data frame.",
        length(brewer_palette_name), " != ", ncol(df),
        collapse = ''
      ))
    }
  } else {
    brewer_palette_name <- rep.int(brewer_palette_name, ncol(df))
  }
  
  metadata_colours <- list()
  extra_colours <- c("darkred","darkblue","black","gray","white")
  for (i in seq(ncol(df))) {
    if (is.null(palette_vals)) {
      # Find unique values
      uniq_vals <- unique(df[,i])
      
      # Determine the number of values to take from the palette
      if (is.null(palette_size))
        draw_from_palette <- max(length(uniq_vals), 3)
      else
        draw_from_palette <- palette_size
      
      var_pal <- c(RColorBrewer::brewer.pal(n = draw_from_palette,
                                            name = brewer_palette_name[i]),
                   extra_colours)
    } else {
      var_pal <- palette_vals
    }
    metadata_colours[i] <- list(setNames(c(var_pal[1:length(uniq_vals)]),
                                         sort(uniq_vals)))
  }
  names(metadata_colours) <- names(df)
  return(metadata_colours)
}



plot_genome_coverage_heatmap <- function(cov_mat,
                                         meta_df,
                                         column_labels = "group",
                                         annotations = NULL,
                                         brewer_palette_name = c("Dark2", "Set3"),
                                         legend_title = "Mean\ncoverage") {
  n_sig_digits <- 4
  hm_legend <- list(
    title = legend_title,
    at = round(c(
      min(cov_mat),
      mean(cov_mat),
      max(cov_mat)),
      n_sig_digits)
  )
  
  if (round(max(cov_mat), n_sig_digits) == 0) {
    stop("Heatmap legend significant digits results in 0")
  }

  if (!is.null(annotations)) {
    annot_palette <- meta_df %>% 
      dplyr::select(tidyselect::all_of(annotations)) %>% 
      colour_metadata(brewer_palette_name = brewer_palette_name)
  } else {
    annot_palette <- colour_metadata(meta_df,
                                     brewer_palette_name = brewer_palette_name)
  }

  col_ha <- ComplexHeatmap::columnAnnotation(
    df = meta_df %>% 
      dplyr::select(tidyselect::all_of(annotations)),
    col = annot_palette
  )
  
  ComplexHeatmap::Heatmap(
    matrix = cov_mat,
    top_annotation = col_ha,
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    border = FALSE,
    col = viridis::viridis(100),
    column_labels = meta_df[[column_labels]],
    heatmap_legend_param = hm_legend,
    clustering_distance_columns = "canberra",
    clustering_distance_rows = "canberra"
  ) %>% 
    return()
}


# row_dist = dist(all_coverage_info_agg_dcast_df, method="canberra")
# row_cluster = hclust(row_dist)
# row_dendro = as.dendrogram(row_cluster)
# row_dendro = rotate(row_dendro,c(1:8,10:17,9))
# col_dist = dist(t(all_coverage_info_agg_dcast_df), method="canberra")
# col_cluster = hclust(col_dist)
# col_dendro = as.dendrogram(col_cluster)
# 
# heatmap3(all_coverage_info_agg_dcast_df, Colv = col_dendro, Rowv = rev(row_dendro), scale="none",margins = c(15, 12),col = colorRampPalette(viridis(6),bias=1.4)(1024))
# 
# # Z-score columns (libraries)
# all_coverage_info_agg_dcast_zscored = scale(all_coverage_info_agg_dcast, center = T, scale = T)
# row_dist = dist(all_coverage_info_agg_dcast_zscored, method="canberra")
# row_cluster = hclust(row_dist)
# row_dendro = as.dendrogram(row_cluster)
# row_dendro = rotate(row_dendro,c(3:12,1,2,16,17,15,13:14))
# col_dist = dist(t(all_coverage_info_agg_dcast_zscored), method="canberra")
# col_cluster = hclust(col_dist)
# col_dendro = as.dendrogram(col_cluster)
# col_dendro = rotate(col_dendro,c(3:18,1,2,19:21))
# 
# # Plot heatmap of z-scored libraries by features
# heatmap3(all_coverage_info_agg_dcast_zscored, Colv = col_dendro, Rowv = rev(row_dendro), scale="none",margins = c(15, 12),col = colorRampPalette(viridis(6),bias=1.4)(1024))
# 
# #library_pvalues = parPvclust(cl=NULL, all_coverage_info_agg_dcast_zscored, method.hclust = "average",method.dist = "canberra", nboot = 1000,iseed = 3333)
# feature_pvalues = parPvclust(cl=NULL, t(all_coverage_info_agg_dcast_zscored), method.hclust = "average",method.dist = "canberra", nboot = 1000,iseed = 3333)
# 
# plot(feature_pvalues, hang = -1, cex = 0.5)

# Remelt
# all_coverage_info_agg_dcast_zscored_remelt = melt(all_coverage_info_agg_dcast_zscored)
# 
# ggplot(all_coverage_info_agg_dcast_zscored_remelt, aes(x=Var2, y=Var1, fill = value))+
#       geom_tile(color = "black")+
#       scale_fill_gradient(low = "gray90", high = "black", na.value="white") +
#       xlab("") +
#       ylab("") +
#       theme_minimal() +
#       coord_fixed(ratio = 1:1) +
#       theme(axis.text.x = element_text(color="black",angle = 60, vjust = 1, size = 9, hjust = 1),
#             axis.text.y = element_text(size=9,color="black"))

  
