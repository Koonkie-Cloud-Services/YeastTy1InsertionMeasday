
p.data_dir <- file.path("KOMB", "data")
p.genome_cov_dir <- file.path(p.data_dir, "feature_heatmapping")
p.input_genome_coverage_file_names <- c(
  "feature-postive"="all_positive_features.tsv",
  "feature-negative"="all_negative_features.tsv",
  "intergenic-positive"="all_positive_antifeatures.tsv",
  "intergenic-negative"="all_negative_antifeatures.tsv"
)

.load_genome_coverage_table <- function(file_path, col_name_lookup) {
  read.table(file = file_path,
             header = FALSE,
             row.names = NULL,
             sep = "\t") %>%
    dplyr::rename(tidyselect::any_of(col_name_lookup)) %>% 
    return()
}

load_genome_coverage_data_frame <- function(file_path, index) {
  col_name_lookup <- c(
    "library"="V1",
    "chr"="V2",
    "start"="V3",
    "end"="V4"
  )
  if ( stringr::str_detect(string = index, pattern = "feature")) {
    df <- .load_genome_coverage_table(
      file_path,
      col_name_lookup = c(col_name_lookup,
                          "feature.name"="V5",
                          "strand"="V6",
                          "coverage"="V7",
                          "feature.length"="V8")
    )
  } else if ( stringr::str_detect(string = index, pattern = "intergenic")) {
    df <- .load_genome_coverage_table(
      file_path,
      col_name_lookup = c(col_name_lookup,
                          "coverage"="V5",
                          "feature.length"="V6")
    ) %>% 
      dplyr::mutate(index = index) %>% 
      tidyr::separate(col = index,
                      into = c("feature.name", "strand")) %>% 
      dplyr::mutate(strand = ifelse(strand == "positive",
                                    "+",
                                    "-"))
  }
  df %>% 
    dplyr::mutate(mean_cov = coverage/feature.length) %>% 
    return()
}

cov_df <- purrr::imap(
  .x = setNames(
    file.path(p.genome_cov_dir, p.input_genome_coverage_file_names),
    names(p.input_genome_coverage_file_names)
  ),
  .f = load_genome_coverage_data_frame) %>% 
  dplyr::bind_rows()
