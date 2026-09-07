# lollipop function

lollipop_plot <- function(data, condition, var_A, var_B, order = TRUE, order_vec, col_A, col_B, zscore_data = FALSE) {
  condition_sym <- rlang::sym(condition)
  volcano_df <- data %>%
    filter(!is.na(!!condition_sym)) %>%
    group_by(name) %>%
    mutate(value_z = if (isTRUE(zscore_data)) as.numeric(scale(value)) else as.numeric(value)) %>%
    summarise(
      delta_z = mean(value_z[!!condition_sym == var_A], na.rm = TRUE) -
        mean(value_z[!!condition_sym == var_B], na.rm = TRUE),
      pvalue  = t.test(value_z[!!condition_sym == var_A],
                       value_z[!!condition_sym == var_B])$p.value,
      .groups = "drop"
    ) %>%
    mutate(
      padj        = p.adjust(pvalue, method = "BH"),
      neg_log10_p = -log10(pvalue),
      direction   = ifelse(delta_z > 0, paste("Up in", var_A), paste("Up in", var_B)),
      sig         = ifelse(padj < 0.1, "FDR < 0.1", ifelse(pvalue < 0.05, "P < 0.05", "Non"))
    ) %>%
    arrange(delta_z)

  if (order == TRUE) {
    volcano_df <- volcano_df %>%
      mutate(name = factor(name, levels = name[order(delta_z)]))
    print(volcano_df$name)
  } else {
    volcano_df <- volcano_df %>%
      mutate(name = factor(name, levels = order_vec))
  }

  dir_A <- paste("Up in", var_A)
  dir_B <- paste("Up in", var_B)

  # ensure the legend shows a break at the true min (and max) of the data,
  # in addition to the automatically chosen "pretty" breaks
  round_up   <- function(x, digits = 1) ceiling(x * 10^digits) / 10^digits
  round_down <- function(x, digits = 1) floor(x * 10^digits) / 10^digits

  size_range  <- range(volcano_df$neg_log10_p, na.rm = TRUE)
  size_breaks <- seq(size_range[1], size_range[2], length.out = 4)

  n <- length(size_breaks)
  size_breaks[1]        <- round_up(size_breaks[1], 1)      # stays >= true min - safe
  size_breaks[n]         <- round_down(size_breaks[n], 1)    # stays <= true max - safe
  size_breaks[-c(1, n)] <- round(size_breaks[-c(1, n)], 1)   # interior points, normal rounding is fine
  size_breaks <- unique(size_breaks)

  p2 <- ggplot(volcano_df, aes(x = delta_z, y = name)) +
    geom_segment(aes(x = 0, xend = delta_z, y = name, yend = name), color = "grey70") +
    geom_point(aes(color = direction, size = neg_log10_p)) +
    scale_color_manual(values = setNames(c(col_A, col_B), c(dir_A, dir_B))) +
    ggnewscale::new_scale_color() +
    geom_point(
      data = subset(volcano_df, sig != "Non"),
      aes(x = delta_z, y = name, size = neg_log10_p, color = sig),
      shape = 21, fill = NA, stroke = 1
    ) +
    scale_color_manual(name = "Significance",
                       values = c("FDR < 0.1" = "black", "P < 0.05" = "grey60")) +
    scale_size_continuous(name = "-log10(p-value)", breaks = size_breaks) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    labs(x = paste0("Δ z-score (", var_A, " - ", var_B, ")"), y = NULL) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 14),
      axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
      axis.text.y  = element_text(size = 14),
      legend.margin  = margin(t = 0),
      legend.spacing = unit(0.2, "cm"),
      legend.text = element_text(size = 12)
    ) +
    ggtitle("Imaging modules")

  return(list(data = volcano_df, plot = p2))
}
