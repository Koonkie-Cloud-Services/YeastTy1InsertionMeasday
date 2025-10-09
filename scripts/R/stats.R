
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
      "wt.mean"="estimate1",
      "mut.mean"="estimate2"
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
        t_test(df=.,
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
