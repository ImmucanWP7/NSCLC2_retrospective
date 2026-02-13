read_vcf <- function(x){
  vcf <- read.vcfR(file = x, verbose = FALSE)
  # if not mutations are found in the vcf file
  if(nrow(vcf@fix) == 0){
    return(NULL)
  }
  vcf <- vcfR2tidy(vcf,format_fields = c("AF", "DP", "AD"))
  # get mapping for this patient

  cur_dat <- vcf$gt
  variants <- vcf$fix

  # split column with annotations by '|'
  ann_split <- str_split(variants$ANN,"\\|")

  # create final data frame
  df_variants <- variants |> dplyr::select(ChromKey, CHROM, POS, REF, ALT) |>
    left_join(cur_dat) |>
    left_join(variants |> dplyr::select(CHROM, POS, REF, ALT, COSMIC, CancerHotspots, gnomAD)) |>
    dplyr::mutate(gene = map_chr(ann_split, c(4)),
                  type = map_chr(ann_split, c(2)),
                  transcript = map_chr(ann_split, c(7)),
                  HGVS_C = map_chr(ann_split, c(10)),
                  HGVS_P = map_chr(ann_split, c(11)))

  # add patient and sample id
  df_variants$immucan_id <- str_split_i(basename(x),"-FIXT",1)
  df_variants$Sample <- str_split_i(basename(x),"_vs_",1)

  return(df_variants)
}
