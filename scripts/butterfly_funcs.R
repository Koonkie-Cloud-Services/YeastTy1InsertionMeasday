
load_trna_region_data_frame <- function(file_path) {
  col_name_lookup = c(
    "trna.name"="name"
  )
  readr::read_table(file = file_path,
                    col_names = TRUE, 
                    show_col_types = FALSE,
                    progress = FALSE) %>%
    dplyr::rename(tidyselect::any_of(col_name_lookup)) %>% 
    dplyr::mutate(
      chr = gsub(pattern = 'chr', replacement = '', x = chr),
      chr = as.roman(chr),
      chr = paste0("chr", chr)
    )
}

.load_coverage_matrix <- function(file_path) {
  col_name_lookup = c(
    "rel.position"="rel_position"
  )
  read.table(file = file_path,
             header = TRUE,
             row.names = NULL,
             sep = "\t") %>%
    dplyr::rename(tidyselect::any_of(col_name_lookup)) %>% 
    return()
}

.add_trna_names <- function(df, trna_feat_df) {
  df %>% 
    group_by(chr, rel.position, proximity, strand) %>%
    mutate(trna = row_number()) %>%
    group_by(chr, trna) %>%
    mutate(
      trna.coord = case_when(
        strand == '-' & proximity == "upstream" ~ min(position),
        strand == '+' & proximity == "upstream" ~ max(position) + 2,
        strand == '-' & proximity == "downstream" ~ max(position),
        strand == '+' & proximity == "downstream" ~ min(position),
        .default = NA)
    ) %>% 
    ungroup() %>%
    left_join(trna_feat_df %>%
                distinct(chr, trna.name, start, end) %>%
                pivot_longer(where(is.numeric),
                             names_to = "position",
                             values_to = "trna.coord") %>% 
                select(!all_of(c("position"))),
              by = c("chr", "trna.coord")) %>% 
    dplyr::relocate(any_of(c("rel.position", "trna.name", "strand")),
                    .after = chr) %>% 
    dplyr::select(!any_of(c("trna", "trna.coord", "start", "end", "score"))) %>%
    return()
}

.separate_strand_and_proximity <- function(df, index) {
  df %>% 
    dplyr::mutate(
      strand_pos = index
    ) %>% 
    tidyr::separate(strand_pos,
                    into = c("proximity", "strand"),
                    sep = '-') %>%
    dplyr::mutate(strand = case_when(
      strand == "positive" ~ '+',
      strand == "negative" ~ '-',
      .default = NA
    )) %>% 
    dplyr::relocate(any_of(c("rel.position", "trna.name", "proximity", "strand")),
                    .after = chr) %>% 
    return()
}

.pivot_coverage_longer_and_aggregate <- function(df) {
  df %>% 
    dplyr::mutate(
      rel.position.char = as.character(rel.position)
    ) %>%
    dplyr::select(!all_of(c("position", "rel.position"))) %>%
    tidyr::pivot_longer(where(is.numeric),
                        names_to = "library",
                        values_to = "coverage") %>% 
    dplyr::mutate(rel.position = as.numeric(rel.position.char)) %>% 
    dplyr::group_by(library, chr, trna.name, proximity, strand, rel.position) %>%
    summarise(
      coverage = sum(coverage),
      .groups = "keep"
    ) %>% 
    ungroup() %>% 
    return()
}

load_coverage_data_frame <- function(file_path, index, ...) {
  .load_coverage_matrix(file_path) %>% 
    .separate_strand_and_proximity(index = index) %>% 
    .add_trna_names(...) %>% 
    .pivot_coverage_longer_and_aggregate() %>%
    return()
}

normalise_sequencing_library_coverage <- function(df, read_counts, read_norm_type) {
  df %>% 
    left_join(read_counts %>% 
                filter(norm_type == read_norm_type),
              by = "library") %>%
    mutate(coverage.norm = coverage/seq_library_size) %>% 
    select(!all_of(c("seq_library_size", "norm_type"))) %>% 
    return()
}


.get_trna_gene_df <- function(trna_feat_df, start_pos=1000) {
  trna_feat_df %>% 
    dplyr::mutate(length = end - start) %>%
    dplyr::reframe(length = median(length)) %>% 
    dplyr::mutate(xmin = start_pos,
                  xmax = start_pos + length,
                  y = 0,
                  forward = TRUE) %>% 
    return()
}

.get_trna_tss_df <- function(start_pos=1000) {
  data.frame(x = start_pos,
             y = 0,
             label = "TSS",
             forward = TRUE)
}

plot_butterfly <- function(df,
                           trna_feat_df,
                           annot_df,
                           span = 0.75,
                           level = 0.95,
                           colour_by = library,
                           colour_palette = NULL,
                           nucleosome_shape = "circle",
                           .annot_width = 70,
                           .annot_height = NULL,
                           .annot_size = 3,
                           .smooth_se = FALSE,
                           .linewidth = 0.5) {
  
  plt <- df %>%
    mutate(
      rel.position = ifelse(strand == "-",
                            (rel.position - 1000) * -1,
                            rel.position),
      rel.position = ifelse(proximity == "downstream",
                            rel.position + 1000,
                            rel.position),
      coverage.norm = ifelse(strand == "-",
                             coverage.norm * -1,
                             coverage.norm)) %>% 
    ggplot2::ggplot()
  
  if (nucleosome_shape == "rectangle") {
    # Automatically determine the optimal height
    if (is.null(.annot_height)) {
      .annot_height <- max(df$coverage.norm) / 100
    }
    
    plt <- plt +
      ggplot2::geom_rect(
        ggplot2::aes(xmin = x-.annot_width,
                     xmax = x+.annot_width,
                     ymin = y-.annot_height,
                     ymax = y+.annot_height),
        alpha = 0.5,
        colour = "darkgrey",
        fill = "lightgrey",
        data = annot_df
      )
  } else if (nucleosome_shape == "circle") {
    plt <- plt +
      ggplot2::geom_point(
        data = annot_df,
        ggplot2::aes(x=x),
        size = .annot_size,
        shape = 21,
        y = 0,
        colour = "darkgrey",
        fill = "lightgrey"
      )
  } else
    stop(paste("Unexpected nucleosome shape:", nucleosome_shape))
  plt <- plt +
    ggplot2::geom_smooth(
      ggplot2::aes(x = rel.position,
                   y = coverage.norm,
                   group = interaction(library, proximity, strand)),
      colour = "black",
      method = "loess",
      formula = y ~ x,
      span = span,
      se = .smooth_se,
      level = level,
      linewidth = .linewidth * 2) + 
    ggplot2::geom_smooth(
      ggplot2::aes(x = rel.position,
                   y = coverage.norm,
                   group = interaction(library, proximity, strand),
                   colour = {{colour_by}}),
      method = "loess",
      formula = y ~ x,
      span = span,
      se = .smooth_se,
      level = level,
      linewidth = .linewidth) +
    xlab("Relative position from TSS (bp)") +
    ylab("Normalized 5'-most nucleotide coverage") +
    ggplot2::facet_wrap(~subunit, scales = "fixed") +
    scale_x_continuous(breaks = seq(0, 2000, by = 250),
                       labels = seq(-1000, 1000, by=250)) +
    gggenes::geom_feature(ggplot2::aes(x = x,
                                       y = y,
                                       forward = forward),
                          data = .get_trna_tss_df()) +
    gggenes::geom_feature_label(ggplot2::aes(x = x,
                                             y = y,
                                             forward = forward,
                                             label = label),
                                data = .get_trna_tss_df()) +
    gggenes::geom_gene_arrow(data = .get_trna_gene_df(trna_feat_df),
                             ggplot2::aes(xmin = xmin,
                                          xmax = xmax,
                                          forward = forward,
                                          y = y)) +
    guides(color = guide_legend(title = "Library",
                                override.aes = list(linewidth = 2),
                                ncol = 1)) +
    theme(
      axis.text.y = element_text(color= "black", size=10),
      axis.text.x = element_text(color= "black", size=10, angle=45, hjust = 1),
      legend.position = "right",
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_line(colour="gray", linewidth = 0.1),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "gray")
    )
  
  if (!is.null(colour_palette)) {
    plt <- plt +
      scale_colour_manual(values = colour_palette)
  }
  
  return(plt)
}
