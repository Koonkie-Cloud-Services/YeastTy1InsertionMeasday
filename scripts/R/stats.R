
t_test <- function(df, values_col, ...) {
  t.test(
    x = df %>% 
      dplyr::filter(strain.type == "Wild-type") %>% 
      dplyr::pull({{values_col}}), 
    y = df %>% 
      dplyr::filter(strain.type == "Mutant") %>% 
      dplyr::pull({{values_col}}),
    ...
  ) %>% 
    broom::tidy() %>% 
    dplyr::rename(
      "wt.mean" = "estimate1",
      "mut.mean" = "estimate2"
    ) %>% 
    return()
}

subunit_wise_t_test_helper <- function(df, ...) {
  # Required columns: 
  # * strain.type: either Mutant or Wild-type
  # * group: the grouping variable to perform comparisons between MT and WT
  unique_groups <- df %>%
    filter(strain.type != "Wild-type") %>% 
    pull(group) %>% 
    unique()
  
  res_cat_df <- data.frame()
  for (group_name in unique_groups) {
    subgroup_df <- df %>% 
      filter(group == group_name | group == "WT")
    obs_matrix <- subgroup_df %>% 
      group_by(strain.type) %>% 
      tally()
    
    if (min(obs_matrix$n) < 3) {
      message(paste0("Skipping t-test for '", group_name, "' due to insufficient observations."))
    } else if ( nrow(obs_matrix) != 2 ) {
      message(paste0("Both mutant and wild-type were not found for '", group_name, "'."))
    } else {
      res_df <- subgroup_df %>% 
        t_test(df = .,
               ...) %>% 
        dplyr::mutate(group = group_name)
      if (nrow(res_cat_df) == 0)
        res_cat_df <- res_df
      else
        res_cat_df <- bind_rows(res_cat_df, res_df)
    }
  }
  return(res_cat_df)
}

sugar_wise_t_test_helper <- function(df, values_col, ...) {
  # Required columns:
  # * sugar: either "dextrose" or "Galactose"
  # * group: the grouping variable to perform comparisons between dextrose and WT Galactose
  unique_groups <- df %>%
    filter(sugar == "Dextrose") %>%
    pull(group) %>%
    unique()

  wt_Galactose_df <- df %>%
    filter(strain.type == "Wild-type", sugar == "Galactose")

  res_cat_df <- data.frame()
  for (group_name in unique_groups) {
    subgroup_df <- bind_rows(
      df %>% filter(group == group_name, sugar == "Dextrose"),
      wt_Galactose_df
    )

    n_dextrose <- subgroup_df %>% filter(sugar == "Dextrose") %>% nrow()
    n_Galactose <- subgroup_df %>% filter(sugar == "Galactose") %>% nrow()

    if (min(n_dextrose, n_Galactose) < 3) {
      message(paste0("Skipping t-test for '", group_name, "' due to insufficient observations."))
    } else {
      res_df <- t.test(
        x = wt_Galactose_df %>% dplyr::pull({{values_col}}),
        y = subgroup_df %>% filter(sugar == "Dextrose") %>% dplyr::pull({{values_col}}),
        ...
      ) %>%
        broom::tidy() %>%
        dplyr::rename(
          "wt.galactose.mean" = "estimate1",
          "dextrose.mean" = "estimate2"
        ) %>%
        dplyr::mutate(group = group_name)

      if (nrow(res_cat_df) == 0)
        res_cat_df <- res_df
      else
        res_cat_df <- bind_rows(res_cat_df, res_df)
    }
  }
  return(res_cat_df)
}


paired_log_cpm_wilcoxon <- function(df, ref_library, test_library,
                                    count_col = "coverage.count",
                                    size_col = "seq_library_size",
                                    pseudocount = 0.5,
                                    conf.level = 0.99,
                                    id_col = "trna.name",
                                    exact = FALSE) {
  # Paired (by tRNA gene) Wilcoxon signed-rank test on log2 counts-per-million.
  # Required columns: library, <id_col> (unique per library), <count_col> (raw counts),
  # <size_col> (library size).
  # Counts are over-dispersed (variance >> mean) and contain a few zeros, so
  # log2((count + pseudocount) / library size * 1e6) stabilises the variance, keeps zeros
  # finite, and makes the paired difference a log2 fold change. The rank-based test does not
  # assume normality of the differences. Genes with no coverage in either library carry no
  # information and are removed.
  to_cpm <- function(lib) {
    df %>%
      dplyr::filter(library == lib) %>%
      dplyr::transmute(
        id = .data[[id_col]],
        count = .data[[count_col]],
        cpm = .data[[count_col]] / .data[[size_col]] * 1e6,
        log2.cpm = log2((.data[[count_col]] + pseudocount) / .data[[size_col]] * 1e6)
      )
  }
  paired_df <- dplyr::inner_join(
    to_cpm(ref_library), to_cpm(test_library),
    by = "id", suffix = c(".ref", ".test")
  ) %>%
    dplyr::filter(count.ref + count.test > 0)

  if (all(paired_df$log2.cpm.test == paired_df$log2.cpm.ref)) {
    return(tibble::tibble(
      n.genes = nrow(paired_df), statistic = NA_real_, p.value = NA_real_,
      log2.fc = NA_real_, conf.low = NA_real_, conf.high = NA_real_,
      wt.mean = mean(paired_df$cpm.ref), mut.mean = mean(paired_df$cpm.test)
    ))
  }

  wt <- suppressWarnings(
    stats::wilcox.test(
      x = paired_df$log2.cpm.test,
      y = paired_df$log2.cpm.ref,
      paired = TRUE,
      conf.int = TRUE,
      conf.level = conf.level,
      exact = exact
    )
  )
  tibble::tibble(
    n.genes = nrow(paired_df),
    statistic = unname(wt$statistic),
    p.value = wt$p.value,
    log2.fc = unname(wt$estimate),  # Hodges-Lehmann median paired log2 fold change
    conf.low = wt$conf.int[1],
    conf.high = wt$conf.int[2],
    wt.mean = mean(paired_df$cpm.ref),
    mut.mean = mean(paired_df$cpm.test)
  )
}
