IFInteractionCounts <- function(data_path, graph_type = "expansion",out_folder, n_workers = 2) {

  cur_files <- list.files(data_path,
                          pattern = ".rds$",
                          full.names = TRUE)

  bplapply(cur_files,function(x){
    cur_sce <- readRDS(x)
    print(unique(cur_sce$sample_id))
    # print(cur_sce)

    if (ncol(cur_sce) > 20000) {

      tumor_sce <- cur_sce[,cur_sce$dist_to_tumor < 500]

      # the average cell radius is roughly 6.6-6.8 µm. we will therefore use 17 µm as a distance for graph building to have roughly 4 µm between cells
      tumor_sce <- buildSpatialGraph(tumor_sce, img_id = "sample_id",
                                     type = graph_type, threshold = 17, coords = c("nucleus.x", "nucleus.y"))

      count_df <- countInteractions(tumor_sce, group_by = "sample_id",
                                    label = "celltype",
                                    colPairName = "expansion_interaction_graph")

      # compute total number of cells
      ct_count <- table(tumor_sce$celltype)

      count_df <- count_df %>%
        as.data.frame() %>%
        mutate(abs_count = ct_count[as.character(from_label)],
               total_interactions = ct * abs_count)

      # save interaction count dataframe
      if(!dir.exists(out_folder)){
        dir.create(out_folder)
      }
      image_id <- unique(cur_sce$sample_id)
      write.csv(count_df,paste0(out_folder,image_id,".csv"),row.names = FALSE)

    }
  }, BPPARAM = MulticoreParam(workers = n_workers))
}
