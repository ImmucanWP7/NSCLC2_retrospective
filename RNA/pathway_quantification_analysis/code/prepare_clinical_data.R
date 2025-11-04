
# select the samples for which RNA is available
ind_av_samples = which(clinical_data$immucan_id %in% colnames(NSCLC2))
subset_clinical_data = clinical_data[ind_av_samples,]

# rename LCNEC
ind_LCNEC <- which(subset_clinical_data$simple_histology == "Large cell neuroendocrine carcinoma (LCNEC)")
subset_clinical_data$simple_histology[ind_LCNEC] <- subset_clinical_data$simple_histology[ind_LCNEC] <- "LCNEC"

# keep the most important variables
subset_clinical_data <- subset_clinical_data[,c(2, 7, 27, 29)]
rownames(subset_clinical_data) <- subset_clinical_data$immucan_id
subset_clinical_data$immucan_id <- NULL

# specify the colors for the clinical variables
annotation_colors = list(
  "gender" = c("Male"="black", "Female"="white"),
  "simple_histology" = c("Adenocarcinoma"="blue", "Adenosquamous carcinoma"="purple", "LCNEC"="black", "Other"="white", "Squamous cell carcinoma"="red", "unknown"="grey"),
  "simple_stage" = c("I"="white", "II"="gray", "III"= "brown", "IV"="black")
)
