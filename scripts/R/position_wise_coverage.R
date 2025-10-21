buffer_length_bp <- 500

prep_position_wise_feature_coverage_df <- function(df,
                                                   gene_df,
                                                   min_gene_read_count = 1,
                                                   min_gene_coverage = 0,
                                                   min_feature_length = buffer_length_bp * 2,
                                                   add_zero_coverage_positions = TRUE,
                                                   summarise_by_position = TRUE) {
  input_gene_count <- nrow(gene_df)
  feature_length_list <- gene_df %>% 
    filter(length >= min_feature_length) %>% 
    pull(feature)
  message(
    paste("Dropped", input_gene_count - length(feature_length_list), "genes shorter than", min_feature_length, "bp")
  )
  
  feature_cov_list <- df %>% 
    group_by(library, feature) %>%
    dplyr::summarise(
      read_total = sum(coverage.raw),
      coverage_total = sum(coverage),
      .groups = "keep"
    ) %>% 
    filter(
      read_total >= min_gene_read_count,
      coverage_total >= min_gene_coverage
    ) %>% 
    pull(feature) %>% 
    unique()
  message(
    paste("Dropped", input_gene_count - length(feature_cov_list), "low-coverage genes")
  )
  
  feature_list <- base::intersect(feature_length_list, feature_cov_list)
  message(
    paste("Analyzing ", length(feature_list), "/",  input_gene_count, " genes", sep = "")
  )
  
  result_df <- df %>%
    filter(
      feature %in% feature_list
    ) %>%
    left_join(gene_df,
              by = "feature") %>% 
    mutate(
      rel_position = case_when(
        start + buffer_length_bp > position ~ position - start,
        stop - buffer_length_bp < position ~ position - stop + (2 * buffer_length_bp),
        .default = NA
      ),
      rel_position = case_when(
        strand == '-' ~ (rel_position * -1) +  (2 * buffer_length_bp),
        .default = rel_position
      )
    ) %>%
    dplyr::filter(
      !is.na(rel_position)
    ) %>%
    select(library, seqid, rel_position, strand, feature, coverage)
  
  if (add_zero_coverage_positions == TRUE) {
    result_df <- result_df %>%
      dplyr::group_by(library, seqid, strand, feature) %>%
      tidyr::expand(rel_position = seq.int(from = buffer_length_bp * -1,
                                           to = buffer_length_bp * 3)) %>% 
      dplyr::mutate(coverage = NA) %>% 
      left_join(result_df,
                by = c("library", "seqid", "strand", "feature", "rel_position"),
                suffix = c("_df1", "_df2")) %>%
      mutate(
        coverage = coalesce(coverage_df1, coverage_df2),
        coverage = coalesce(coverage, 0)
      ) %>%
      select(!all_of(c("coverage_df1", "coverage_df2")))
  }
  
  if (summarise_by_position == TRUE) {
    result_df <- result_df %>% 
      group_by(library, rel_position) %>% 
      dplyr::reframe(
        coverage.sum = sum(coverage),
        coverage.mean = mean(coverage),
        coverage.median = median(coverage)
      )
  }
  
  if (nrow(result_df) == 0)
    stop("All genes dropped by filtering thresholds.")
  return(result_df)
}

plot_position_coverage <- function(df,
                                   palette_vals,
                                   plot_type = "smooth",
                                   label_y_position = 0,
                                   text_gap_size = 100,
                                   loess_span = 0.01,
                                   x_gap_size = 20,
                                   x_scale = "neighbourhood",
                                   free_y_scale = FALSE,
                                   y_max = NULL) {
  if ( ! plot_type %in% c("smooth", "line", "bar", "density") )
    stop("Invalid plot type")
  
  plt <- df %>% 
    ggplot(aes(
      x = rel_position,
      y = coverage,
      group = library
    )) +
    ggplot2::scale_colour_manual(values = palette_vals) +
    ggplot2::guides(colour = guide_legend(title = "")) +
    theme_bw() +
    theme(plot.margin = unit(c(2,1,1,1), "lines")) +
    ggplot2::ylab("Relative number of Ty1 insertions\n(normalised coverage)")
  
  if (x_scale == "percent") {
    plt <- plt +
      ggplot2::geom_vline(xintercept = 0,
                          linetype = "dashed") +
      ggplot2::geom_vline(xintercept = 100,
                          linetype = "dashed") +
      ggplot2::scale_x_continuous(
        breaks = c(buffer_length_bp * -1, 0, 100, buffer_length_bp + 100),
        labels = c(buffer_length_bp * -1, 0, 100, buffer_length_bp)
      ) +
      ggplot2::xlab("Percentage of way through ORF")
  } else {
    plt <- plt +
      ggplot2::geom_vline(aes(xintercept = min(rel_position) + buffer_length_bp),
                          linetype = "dashed") +
      ggplot2::geom_vline(aes(xintercept = max(rel_position) - buffer_length_bp),
                          linetype = "dashed") +
      ggplot2::scale_x_continuous(
        breaks = c(buffer_length_bp * -1, 0, buffer_length_bp, buffer_length_bp * 2, buffer_length_bp * 3),
        labels = c(buffer_length_bp * -1, 0, buffer_length_bp, 0, buffer_length_bp)
      )
  }
  
  if (free_y_scale == TRUE) {
    plt <- plt +
      ggplot2::facet_grid(rows = vars(subunit),
                          scales = "free")
  } else {
    plt <- plt +
      ggplot2::facet_grid(rows = vars(subunit)) +
      ggplot2::annotate("text",
                        x = 0 - text_gap_size,
                        y = label_y_position,
                        size = 3,
                        label = "TSS") +
      ggplot2::annotate("text",
                        x = 1000 + text_gap_size,
                        y = label_y_position,
                        size = 3,
                        label = "TTS") +
      ggplot2::coord_cartesian(clip = 'on',
                               ylim = c(0, 7E-6)) +
      ggplot2::scale_y_continuous(
        breaks = c(0, 7E-6)
      )
  }
  
  if (plot_type == "smooth") {
    plt <- plt +
      geom_smooth(se = FALSE,
                  colour = "black",
                  method = "loess",
                  formula = y ~ x,
                  span = loess_span,
                  level = 0.98,
                  linewidth = 2) +
      geom_smooth(se = FALSE,
                  method = "loess",
                  formula = y ~ x,
                  span = loess_span,
                  level = 0.98,
                  aes(colour = library),
                  linewidth = 1.5)
  }
  else if (plot_type == "line") {
    plt <- plt +
      geom_line(colour = "black",
                linewidth = 2) +
      geom_line(aes(colour = library),
                linewidth = 1.5)
  }
  else if (plot_type == "bar") {
    plt <- plt +
      geom_bar(
        colour = "black",
        aes(fill = library),
        stat = "identity"
      )
  }
  else {
    stop("Unimplemented plot type.")
  }
  if (!is.null(y_max)) {
    plt <- plt +
      ggplot2::geom_rect(
        xmin = buffer_length_bp - x_gap_size,
        xmax = buffer_length_bp + x_gap_size,
        ymin = 0,
        ymax = y_max,
        fill = "white",
        colour = "white",
        alpha = 0.8,
        inherit.aes = FALSE
      )
  }
  return(plt)
}
