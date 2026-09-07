## helper functions
### simple boxplot
bplot <- function(ct, markers, mutation, rownum, histology, stat_size = 2.5){
  mean_exp %>%
    filter(celltypes == ct,
           name %in% markers) %>%
    left_join(traits, by = "immucan_id") %>%
    filter(simple_histology %in% histology,
           !is.na(.data[[mutation]])) %>%
    ggplot(aes(x = .data[[mutation]], y = value,fill= .data[[mutation]]))+
    geom_boxplot(outliers = FALSE)+
    geom_jitter(aes(color = .data[[mutation]]),position = position_jitterdodge(jitter.width = 0.2), size= 0.5)+
    scale_fill_manual(values = alpha(c("red", "blue"),alpha = 0.6))+
    scale_color_manual(values = c("red", "blue"))+
    expand_limits(y = 0) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    facet_wrap(.~name, scales= "free", nrow = rownum)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
    stat_compare_means(size = stat_size, label = "p.format")+
    theme_bw()+
    theme(axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 6, angle = 90, vjust = 1, hjust = 0.5),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6),
          strip.text = element_text(size = 6),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7))+
    theme(legend.position = "none")
}


### boxplot with p value correction


bplot_cor <- function(ct,
                      markers,
                      mutation,
                      histology  = "Adenocarcinoma",
                      stat_size  = 2,
                      rownum     = 1,
                      sig_format = FALSE,
                      p_adjust   = "BH",
                      fdr_cut    = 0.1,      # only label comparisons passing this
                      digits     = 3,
                      headroom   = 0.1) {

  df <- mean_exp %>%
    filter(celltypes == ct, name %in% markers) %>%
    left_join(traits, by = "immucan_id") %>%
    filter(simple_histology %in% histology, !is.na(.data[[mutation]])) %>%
    mutate(across(all_of(mutation), ~ droplevels(factor(.x))))

  stopifnot(nrow(df) > 0)
  if (dplyr::n_distinct(df[[mutation]]) < 2)
    stop("`", mutation, "` has fewer than 2 levels after filtering.")

  stat.test <- df %>%
    group_by(name) %>%
    wilcox_test(reformulate(mutation, response = "value")) %>%
    ungroup() %>%
    adjust_pvalue(method = p_adjust) %>%
    add_significance(
      p.col      = "p.adj",
      output.col = "p.adj.signif",
      cutpoints  = c(0, 0.001, 0.01, 0.05, 0.1, 1),
      symbols    = c("****", "***", "**", "*", "ns")
    )

  # label text: rounded raw p, with a floor so small values don't print as "0.00"
  floor_p <- 10^(-digits)
  stat.test <- stat.test %>%
    mutate(p.label = ifelse(
      p < floor_p,
      sprintf("p < %s", format(floor_p, scientific = FALSE)),
      sprintf("p = %.*f", digits, p)
    ))

  # keep only comparisons that survive FDR
  stat.sig <- stat.test %>% filter(!is.na(p.adj), p.adj < fdr_cut)

  # place labels relative to each panel's own maximum
  if (nrow(stat.sig) > 0) {
    ymax <- df %>%
      group_by(name) %>%
      summarise(ymax = max(value, na.rm = TRUE), .groups = "drop")

    stat.sig <- stat.sig %>%
      left_join(ymax, by = "name") %>%
      group_by(name) %>%
      mutate(y.position = ymax * (1 + headroom * seq_len(n()))) %>%
      ungroup()
  }

  lab <- if (isTRUE(sig_format)) "{p.adj.signif}" else "{p.label}"

  p <- ggplot(df, aes(x = .data[[mutation]], y = value, fill = .data[[mutation]])) +
    geom_boxplot(outliers = FALSE) +
    geom_jitter(aes(color = .data[[mutation]]),
                position = position_jitterdodge(jitter.width = 0.2), size = 0.5) +
    scale_fill_manual(values = alpha(c("red", "blue", "darkgreen"), alpha = 0.6)) +
    scale_color_manual(values = c("red", "blue", "darkgreen")) +
    expand_limits(y = 0) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    facet_wrap(. ~ name, scales = "free", nrow = rownum) +
    labs(title = ct, y = "Mean expression") +
    theme_bw() +
    theme(plot.title      = element_text(size = 8),
          axis.text.y     = element_text(size = 6),
          axis.text.x     = element_text(size = 6, angle = 90, vjust = 1, hjust = 0.5),
          strip.text      = element_text(size = 6),
          axis.title.x    = element_text(size = 7),
          axis.title.y    = element_text(size = 7),
          legend.position = "none")

  if (nrow(stat.sig) > 0) {
    p <- p + stat_pvalue_manual(stat.sig,
                                label      = lab,
                                size       = stat_size,
                                tip.length = 0.01)
  } else {
    message("No comparison passed FDR < ", fdr_cut, " — plotting without labels.")
  }

  attr(p, "stat.test") <- stat.test   # full table, including non-significant rows
  p
}
