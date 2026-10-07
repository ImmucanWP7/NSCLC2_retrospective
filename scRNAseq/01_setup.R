DATA_DIR   <- "/mnt/immucan_volume/processed_data/NSCLC2_cohort/data/scRNAseq/data/"
H5AD_FILE  <- file.path(DATA_DIR, "f678fb47-e51b-4dc5-b23f-f9df43a67ee5.h5ad")

# Bulk WES / clinical metadata from the main NSCLC2 project
SCORES_DIR <- "/mnt/immucan_volume/processed_data/NSCLC2_cohort/Figures/Analysis/scores_for_sharing"
ALL_TRAITS <- file.path(SCORES_DIR, "all_traits.csv")

# All output goes here
OUT_DIR <- file.path("/mnt/immucan_volume/processed_data/NSCLC2_cohort/data/scRNAseq/Rout/")

FILTERED_SCE <- file.path(OUT_DIR, "sce_qc_filtered.rds")

# Output of sc_02: EGFR-subset, gene-symbol rownames, fully in-memory
EGFR_SCE <- file.path(OUT_DIR, "sce_egfr_subset.rds")

library(BiocParallel)
N_WORKERS <- max(1L, parallel::detectCores() - 2L)
register(MulticoreParam(N_WORKERS))

library(zellkonverter)
library(SingleCellExperiment)
library(HDF5Array)
library(rhdf5)
library(dplyr)


# read sce data


sce <- zellkonverter::readH5AD(
  "/mnt/immucan_volume/processed_data/NSCLC2_cohort/data/scRNAseq/data/f678fb47-e51b-4dc5-b23f-f9df43a67ee5.h5ad",
  reader   = "R",
  use_hdf5 = TRUE
)


# ── 2. Filter relevant studies from atlas ─────────────────────────────────────
## EGFR data sets:
# He_Fan -> only EGFR mutated
# Wo_Zhou -> do not provide patient metadata -> cannot identify EGFR mutated patients
# Kim_Lee_2020 -> metadata available
# Maynard_Bivona -> metadata available

lee <- readxl::read_xlsx("/mnt/immucan_volume/processed_data/NSCLC2_cohort/data/scRNAseq/data/Kim_Lee_2020_metadata.xlsx", skip = 1)
lee_pats_mutated <- lee %>% 
  filter(EGFR != "na", 
         !`Tissue origins` %in%  c("nLN","nLung","PE")) %>%
  mutate(EGFR_mutation = ifelse(EGFR == "WT", "not mutated", "mutated"),
         `Patient id` = paste0("Kim_Lee_2020_",`Patient id`)) %>%
  select(`Patient id`, EGFR_mutation)


may <- readxl::read_xlsx("/mnt/immucan_volume/processed_data/NSCLC2_cohort/data/scRNAseq/data/Maynard_et_al_patient_metadata.xlsx", sheet = 2)
may_pats <- may %>%
  filter(Histolgy == "Adenocarcinoma",
         grepl("EGFR", `Oncogenic Driver Mutation`),
         `Biopsy Site` == "Lung",
         `Treatment Hx` == "tx naive") %>%
  mutate(EGFR_mutation = ifelse(`Oncogenic Driver Mutation` == "WT", "not mutated", "mutated"),
         `Patient ID` = paste0("Maynard_Bivona_2020_",`Patient ID`)) %>%
  select(`Patient ID`, EGFR_mutation)

cur_dat <- as.data.frame(colData(sce))
unique(cur_dat$EGFR_mutation)
pats <- cur_dat %>%
  filter(EGFR_mutation %in% c("mutated", "not mutated")) %>%
  pull(donor_id) %>%
  unique() %>%
  as.character()

pats <- c(pats, may_pats$`Patient ID`, lee_pats_mutated$`Patient id`)

egfr_sce <- sce[,sce$donor_id %in% pats]

message("\n=== Materialising assays into RAM (dgCMatrix) ===")
for (aname in assayNames(egfr_sce)) {
  message(sprintf("  %s ...", aname))
  assay(egfr_sce, aname, withDimnames = FALSE) <-
    as(assay(egfr_sce, aname), "dgCMatrix")
}
message("All assays now in-memory.")
print(object.size(egfr_sce), units = "Mb")

saveRDS(egfr_sce,"/mnt/immucan_volume/processed_data/NSCLC2_cohort/data/scRNAseq/Rout/egfr_sce.rds")
