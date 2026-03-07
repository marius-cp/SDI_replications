
p_to_stars <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.01, "***",
                ifelse(p > 0.01 & p < 0.05, "**",
                       ifelse(p > 0.05 & p < 0.1, "*", "")
                )
         )
  )
}


isoline_intercept <- function(s, unc){
  intercept <- unc - s
  return(list(intercept=intercept, score = c(unc-intercept)))
}



plot_sdi_comparison_matrix <- function(comparison_dataframes) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggh4x)
  
  reshape_long <- function(df, name) {
    df %>%
      tibble::rownames_to_column("model_1") %>%
      pivot_longer(cols = -model_1, names_to = "model_2", values_to = name)
  }
  
  combined_long_se <- reshape_long(comparison_dataframes$s, "s") %>%
    left_join(reshape_long(comparison_dataframes$mcb, "mcb"), by = c("model_1", "model_2")) %>%
    left_join(reshape_long(comparison_dataframes$dsc, "dsc"), by = c("model_1", "model_2")) %>%
    left_join(reshape_long(comparison_dataframes$s_pval, "s_pval"), by = c("model_1", "model_2")) %>%
    left_join(reshape_long(comparison_dataframes$mcb_pval, "mcb_pval"), by = c("model_1", "model_2")) %>%
    left_join(reshape_long(comparison_dataframes$dsc_pval, "dsc_pval"), by = c("model_1", "model_2")) %>%
    left_join(reshape_long(comparison_dataframes$mcb_null, "mcb_null"), by = c("model_1", "model_2")) %>%
    left_join(reshape_long(comparison_dataframes$MZ_classic, "MZ_classic"), by = c("model_1", "model_2")) %>%
    left_join(reshape_long(comparison_dataframes$dsc_null, "dsc_null"), by = c("model_1", "model_2")) %>%
    mutate(
      model_1 = factor(model_1, levels = unique(model_1)),
      model_2 = factor(model_2, levels = levels(model_1))
    ) %>%
    filter(as.numeric(model_1) >= as.numeric(model_2))
  
  create_label_data <- function(df, type = c("off", "diag")) {
    type <- match.arg(type)
    df <- if (type == "off") df %>% filter(model_1 != model_2) else df %>% filter(model_1 == model_2)
    
    # df <- df %>%
    #   rowwise() %>%
    #   mutate(
    #     math_labels = list(if (type == "off") c("Delta*widehat(S)", "Delta*widehat(MCB)", "Delta*widehat(DSC)")
    #                        else c("p[MCB]^0", "p[MZ]^phantom(0)", "p[DSC]^0")),
    #     value_labels = list(if (type == "off") c(
    #       if (!is.na(s)) sprintf("= %6s", formatC(s, format = "f", digits = 2)) %>% paste0(ifelse(s_pval < 0.05, "*", "")) else "",
    #       if (!is.na(mcb)) sprintf("= %6s", formatC(mcb, format = "f", digits = 2)) %>% paste0(ifelse(mcb_pval < 0.05, "*", "")) else "",
    #       if (!is.na(dsc)) sprintf("= %6s", formatC(dsc, format = "f", digits = 2)) %>% paste0(ifelse(dsc_pval < 0.05, "*", "")) else ""
    #     ) else c(
    #       if (!is.na(mcb_null)) sprintf("= %6s", formatC(pmax(0, mcb_null), format = "f", digits = 2)) else "",
    #       if (!is.na(MZ_classic)) sprintf("= %6s", formatC(MZ_classic, format = "f", digits = 2)) else "",
    #       if (!is.na(dsc_null)) sprintf("= %6s", formatC(pmax(0, dsc_null), format = "f", digits = 2)) else ""
    #     )),
    #     y_vals = list(c(1.3, 1.0, 0.7))
    #   ) %>%
    #   ungroup() %>%
    #   select(model_1, model_2, math_labels, value_labels, y_vals)
    
    df <- df %>%
      rowwise() %>%
      mutate(
        # math_labels = list(if (type == "off") c("Delta*widehat(S)", "Delta*widehat(MCB)", "Delta*widehat(DSC)")
        #                    else c("p[MCB]^0", "p[MZ]^phantom(0)", "p[DSC]^0")),
        math_labels = list(if (type == "off") c("Delta*widehat(S)", "Delta*widehat(MCB)", "Delta*widehat(DSC)")
                           else c("p[MCB]^0", " ", "p[DSC]^0")),
        value_labels = list(if (type == "off") c(
          if (!is.na(s))   paste0(sprintf("= %5s", formatC(s,   format = "f", digits = 2)), p_to_stars(s_pval))   else "",
          if (!is.na(mcb)) paste0(sprintf("= %5s", formatC(mcb, format = "f", digits = 2)), p_to_stars(mcb_pval)) else "",
          if (!is.na(dsc)) paste0(sprintf("= %5s", formatC(dsc, format = "f", digits = 2)), p_to_stars(dsc_pval)) else ""
        ) else c(
          if (!is.na(mcb_null))   sprintf("= %5s", formatC(pmax(0, mcb_null), format = "f", digits = 2)) else "",
          if (!is.na(MZ_classic)) sprintf("= %5s", formatC(MZ_classic,       format = "f", digits = 2)) else "",
          if (!is.na(dsc_null))   sprintf("= %5s", formatC(pmax(0, dsc_null), format = "f", digits = 2)) else ""
        )),
        y_vals = list(c(1.3, 1.0, 0.7))
      ) %>%
      ungroup() %>%
      select(model_1, model_2, math_labels, value_labels, y_vals)
    
    df %>%
      unnest(c(math_labels, value_labels, y_vals)) %>%
      rename(y = y_vals) %>%
      mutate(x = 1)
  }
  
  label_df <- bind_rows(
    create_label_data(combined_long_se, "off"),
    create_label_data(combined_long_se, "diag")
  )
  
  
  
  ggplot(label_df, aes(x = x, y = y)) +
    geom_tile(
      data = combined_long_se %>% filter(model_1 == model_2),
      aes(x = 1, y = 1),
      fill = "grey90", color = NA
    ) +
    geom_text(
      aes(label = math_labels),
      parse = TRUE,
      hjust = 1,
      size = 5,
      nudge_x = -0.02
    ) +
    geom_text(
      aes(label = value_labels),
      parse = FALSE,
      hjust = 0,
      size = 5,
      family = "mono",        # ← ensures fixed-width alignment
      nudge_x = 0.02
    ) +
    facet_grid2(model_1 ~ model_2, switch = "both", render_empty = FALSE) +
    theme_minimal(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(size = 12),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    ) +
    labs(x = NULL, y = NULL, title = NULL)
  
  ggplot(label_df, aes(x = x, y = y)) +
    geom_tile(
      data = combined_long_se %>% filter(model_1 == model_2),
      aes(x = 1, y = 1),
      fill = "grey90", color = NA
    ) +
    geom_text(
      aes(label = math_labels),
      parse = TRUE,
      hjust = 1,
      size = 4.5,
      nudge_x = -0.01,
      family = "mono"  # optional, keeps baseline consistent
    ) +
    geom_text(
      aes(label = value_labels),
      parse = FALSE,
      hjust = 0,
      size = 4.5,
      nudge_x = 0.01,
      family = "mono"
    ) +
    facet_grid2(model_1 ~ model_2, switch = "both", render_empty = FALSE) +
    theme_minimal(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(size = 12),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    ) +
    labs(x = NULL, y = NULL, title = NULL)
  
  ggplot(label_df, aes(x = x, y = y)) +
    geom_tile(
      data = combined_long_se %>% filter(model_1 == model_2),
      aes(x = 1, y = 1),
      fill = "grey90", color = NA
    ) +
    geom_text(
      aes(label = math_labels),
      parse = TRUE,
      hjust = 1,
      size = 4.3,
      nudge_x = -.175,
      #family = "mono"
    ) +
    geom_text(
      aes(label = value_labels),
      parse = FALSE,
      hjust = 0,
      size = 4.3,
      nudge_x = -.15,  # ← no more right-push
      family = "mono"
    ) +
    facet_grid2(model_1 ~ model_2, switch = "both", render_empty = FALSE) +
    theme_minimal(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(size = 12),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    ) +
    labs(x = NULL, y = NULL, title = NULL)
  
  
}



plot_sdi_comparison_matrix_q <- function(comparison_dataframes) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggh4x)
  
  reshape_long <- function(df, name) {
    df %>%
      tibble::rownames_to_column("model_1") %>%
      pivot_longer(cols = -model_1, names_to = "model_2", values_to = name)
  }
  
  combined_long_se <- reshape_long(comparison_dataframes$s, "s") %>%
    left_join(reshape_long(comparison_dataframes$mcb, "mcb"), by = c("model_1", "model_2")) %>%
    left_join(reshape_long(comparison_dataframes$dsc, "dsc"), by = c("model_1", "model_2")) %>%
    left_join(reshape_long(comparison_dataframes$s_pval, "s_pval"), by = c("model_1", "model_2")) %>%
    left_join(reshape_long(comparison_dataframes$mcb_pval, "mcb_pval"), by = c("model_1", "model_2")) %>%
    left_join(reshape_long(comparison_dataframes$dsc_pval, "dsc_pval"), by = c("model_1", "model_2")) %>%
    left_join(reshape_long(comparison_dataframes$mcb_null, "mcb_null"), by = c("model_1", "model_2")) %>%
    left_join(reshape_long(comparison_dataframes$MZ_classic, "MZ_classic"), by = c("model_1", "model_2")) %>%
    left_join(reshape_long(comparison_dataframes$dsc_null, "dsc_null"), by = c("model_1", "model_2")) %>%
    mutate(
      model_1 = factor(model_1, levels = unique(model_1)),
      model_2 = factor(model_2, levels = levels(model_1))
    ) %>%
    filter(as.numeric(model_1) >= as.numeric(model_2))
  
  create_label_data <- function(df, type = c("off", "diag")) {
    type <- match.arg(type)
    df <- if (type == "off") df %>% filter(model_1 != model_2) else df %>% filter(model_1 == model_2)
    
    # df <- df %>%
    #   rowwise() %>%
    #   mutate(
    #     math_labels = list(if (type == "off") c("Delta*widehat(S)", "Delta*widehat(MCB)", "Delta*widehat(DSC)")
    #                        else c("p[MCB]^0", " ", "p[DSC]^0")),
    #     value_labels = list(if (type == "off") c(
    #       if (!is.na(s)) sprintf("= %6s", formatC(s, format = "f", digits = 2)) %>% paste0(ifelse(s_pval < 0.05, "*", "")) else "",
    #       if (!is.na(mcb)) sprintf("= %6s", formatC(mcb, format = "f", digits = 2)) %>% paste0(ifelse(mcb_pval < 0.05, "*", "")) else "",
    #       if (!is.na(dsc)) sprintf("= %6s", formatC(dsc, format = "f", digits = 2)) %>% paste0(ifelse(dsc_pval < 0.05, "*", "")) else ""
    #     ) else c(
    #       if (!is.na(mcb_null)) sprintf("= %6s", formatC(pmax(0, mcb_null), format = "f", digits = 2)) else "",
    #       if (!is.na(MZ_classic)) sprintf("= %6s", formatC(MZ_classic, format = "f", digits = 2)) else "",
    #       if (!is.na(dsc_null)) sprintf("= %6s", formatC(pmax(0, dsc_null), format = "f", digits = 2)) else ""
    #     )),
    #     y_vals = list(c(1.3, 1.0, 0.7))
    #   ) %>%
    #   ungroup() %>%
    #   select(model_1, model_2, math_labels, value_labels, y_vals)
    
    df <- df %>%
      rowwise() %>%
      mutate(
        math_labels = list(if (type == "off") c("Delta*widehat(S)", "Delta*widehat(MCB)", "Delta*widehat(DSC)")
                           else c("p[MCB]^0", " ", "p[DSC]^0")),
        value_labels = list(if (type == "off") c(
          if (!is.na(s))   paste0(sprintf("= %5s", formatC(s,   format = "f", digits = 2)), p_to_stars(s_pval))   else "",
          if (!is.na(mcb)) paste0(sprintf("= %5s", formatC(mcb, format = "f", digits = 2)), p_to_stars(mcb_pval)) else "",
          if (!is.na(dsc)) paste0(sprintf("= %5s", formatC(dsc, format = "f", digits = 2)), p_to_stars(dsc_pval)) else ""
        ) else c(
          if (!is.na(mcb_null))   sprintf("= %5s", formatC(pmax(0, mcb_null), format = "f", digits = 2)) else "",
          if (!is.na(MZ_classic)) sprintf("= %5s", formatC(MZ_classic,       format = "f", digits = 2)) else "",
          if (!is.na(dsc_null))   sprintf("= %5s", formatC(pmax(0, dsc_null), format = "f", digits = 2)) else ""
        )),
        y_vals = list(c(1.3, 1.0, 0.7))
      ) %>%
      ungroup() %>%
      select(model_1, model_2, math_labels, value_labels, y_vals)
    
    df %>%
      unnest(c(math_labels, value_labels, y_vals)) %>%
      rename(y = y_vals) %>%
      mutate(x = 1)
  }
  
  label_df <- bind_rows(
    create_label_data(combined_long_se, "off"),
    create_label_data(combined_long_se, "diag")
  )
  
  ggplot(label_df, aes(x = x, y = y)) +
    geom_tile(
      data = combined_long_se %>% filter(model_1 == model_2),
      aes(x = 1, y = 1),
      fill = "grey90", color = NA
    ) +
    geom_text(
      aes(label = math_labels),
      parse = TRUE,
      hjust = 1,
      size = 5,
      nudge_x = -0.02
    ) +
    geom_text(
      aes(label = value_labels),
      parse = FALSE,
      hjust = 0,
      size = 5,
      family = "mono",        # ← ensures fixed-width alignment
      nudge_x = 0.02
    ) +
    facet_grid2(model_1 ~ model_2, switch = "both", render_empty = FALSE) +
    theme_minimal(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(size = 12),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    ) +
    labs(x = NULL, y = NULL, title = NULL)
  
  ggplot(label_df, aes(x = x, y = y)) +
    geom_tile(
      data = combined_long_se %>% filter(model_1 == model_2),
      aes(x = 1, y = 1),
      fill = "grey90", color = NA
    ) +
    geom_text(
      aes(label = math_labels),
      parse = TRUE,
      hjust = 1,
      size = 4.5,
      nudge_x = -0.01,
      family = "mono"  # optional, keeps baseline consistent
    ) +
    geom_text(
      aes(label = value_labels),
      parse = FALSE,
      hjust = 0,
      size = 4.5,
      nudge_x = 0.01,
      family = "mono"
    ) +
    facet_grid2(model_1 ~ model_2, switch = "both", render_empty = FALSE) +
    theme_minimal(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(size = 12),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    ) +
    labs(x = NULL, y = NULL, title = NULL)
  
  ggplot(label_df, aes(x = x, y = y)) +
    geom_tile(
      data = combined_long_se %>% filter(model_1 == model_2),
      aes(x = 1, y = 1),
      fill = "grey90", color = NA
    ) +
    geom_text(
      aes(label = math_labels),
      parse = TRUE,
      hjust = 1,
      size = 4.3,
      nudge_x = -.175,
      #family = "mono"
    ) +
    geom_text(
      aes(label = value_labels),
      parse = FALSE,
      hjust = 0,
      size = 4.3,
      nudge_x = -.15,  # ← no more right-push
      family = "mono"
    ) +
    facet_grid2(model_1 ~ model_2, switch = "both", render_empty = FALSE) +
    theme_minimal(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(size = 12),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    ) +
    labs(x = NULL, y = NULL, title = NULL)
  
  
}

