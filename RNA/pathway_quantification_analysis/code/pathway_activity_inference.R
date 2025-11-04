
for (i in 1:length(method_pathway))
{
  current_method <- method_pathway[i]
  
  # Run pathway quantification by current_method
  if (current_method == "mlm") {
    res <- run_mlm(mat = NSCLC2, .source = "source", .target = "target", .mor = weight, net = net_progeny, minsize = 5)
    
  } else if (current_method == "wmean") {
    res <- run_wmean(mat = NSCLC2, net = net_progeny, .source = "source", .target = "target", .mor = weight, times = 1000, minsize = 5, seed = 42)
    res <- res %>% dplyr::filter(statistic == "corr_wmean")
    
  } else if (current_method == "ulm") {
    res <- run_ulm(mat = NSCLC2, net = net_progeny, .source = "source", .target = "target", .mor = "weight", minsize = 5)
    
  } else if (current_method == "aucell") {
    res <- run_aucell(mat = NSCLC2, net = net_progeny, .source = "source", .target = "target", minsize = 5, seed = 42)
    
  } else if (current_method == "gsva") {
    res <- run_gsva(mat = NSCLC2, net = net_progeny, .source = "source", .target = "target", kcdf = "Gaussian", minsize = 5)
    
  } else if (current_method == "ora") {
    res <- run_ora(mat = NSCLC2, net = net_progeny, .source = "source", .target = "target", minsize = 5, seed = 42)
    
  } else if (current_method == "viper") {
    res <- run_viper(mat = NSCLC2, net = net_progeny, .source = "source", .target = "target", .mor = "weight", pleiotropy = TRUE, minsize = 5)
    
  } else if (current_method == "wsum") {
    res <- run_wsum(mat = NSCLC2, net = net_progeny, .source = "source", .target = "target", .mor = "weight", times = 1000, minsize = 5, seed = 42)
    res <- res %>% dplyr::filter(statistic == "corr_wsum")
  } 
  
  
  # Transform the output to a wide matrix: samples * pathways
  wide_res <- res %>%
    select(condition, source, score) %>%      
    pivot_wider(
      names_from  = source,                  
      values_from = score                 
    )
  
  mat <- wide_res %>% column_to_rownames("condition") %>% as.matrix() # set row names to condition
    
  # Scale per sample
  res_mat <- scale(mat)
  
  # Choose color palette
  palette_length = 20
  my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)
  
  my_breaks <- c(seq(floor(min(res_mat)), 0, length.out=ceiling(palette_length/2) + 1),
                 seq(0.05, ceiling(max(res_mat)), length.out=floor(palette_length/2)))
  
  # Representation of the quantifications as a pheatmap
  pdf(paste0( GraphOutDir, "/", outputDirPathwayActivity,"/", "pathway_quantification_NSCLC2_", current_method,".pdf"), height=16.5,width=10)
  pheatmap(res_mat, border_color = NA, color=my_color, breaks = my_breaks, fontsize=14, fontsize_row = 10, 
           clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", clustering_method = "average",
           annotation_row=subset_clinical_data, show_rownames=FALSE, annotation_colors = annotation_colors)
  dev.off()
  
  # Write results of the quantifications
  write.table(res_mat, paste0(outputDir,"/", outputDirPathwayActivity,"/PathwayActivity_NSCLC2_", current_method))

}
