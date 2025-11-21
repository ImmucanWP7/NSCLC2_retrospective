# Publication

This repo contains the code for the analysis of the data for the publication from Schulz, Morfouace et al. entitled "Comprehensive characterization of early-stage Non-Small Cell Lung cancers with multi-modal data integration ".
A preprint of the publication can be found on MedRxiv: [Preprint](https://www.medrxiv.org/content/10.1101/2025.11.10.25339890v1)

# Installation and requirements
The code uses R and was run under R 4.2.2. with BioConductor version 3.16. For R installation please follow instructions [here](https://cran.r-project.org/) and for Bioconductor [here](https://bioconductor.org/).
Numerous libraries need to be installed. They can be found at the beginning of each script.

# Data

The data is currently in the process of being deposited and will be available as follows: 
IMC via [IDR](https://idr.openmicroscopy.org/)
mIF via [IDR](https://idr.openmicroscopy.org/)
[WES](https://ega-archive.org/)
[RNAseq](https://ega-archive.org/)
clinical data: will be provided upon request from EORTC, the owner of the clinical data according to [sharing policies](https://www.eortc.org/app/uploads/2023/06/L-01-POL-01.pdf). Request can be submitted [here](https://www.eortc.be/services/forms/erp/request.aspx) 

# Usage
Typically each subfolder from this repository contains numebered analysis scripts. These should be run in the order or numbering as they depend on each other. Additionally, in most scripts a variable called "mount_path" is set, which defines a parent directory under which the data should be stored. This variable can be modified locally (set to the respective path of the downloaded data) in order to run the script.

# Folders

## IMC

1. read IMC panel 1 data
- The data is provided in folders for each patient
3. perform cell typing

## IMC_p2

1. read IMC panel 1 data
- The data is provided in folders for each patient
3. perform cell typing

## RNA

1. prepare metadata
2. RNA count processing
- for each patient one table is provided.
3. QC plots
4. deconvolution or RNAseq data - not included in the paper

## WES

- immucan_nsclc_wes_tableReport.Rmd mutation identification and Oncoplot
  - starts from table reports generated from vcf files in IMMUcan
- mutational_signatures.R generation of mutational signatures

## MOFA

1. data preparation
2. model training

## Figures

- should be run sequentially in the order of the numbered scripts as script 01 and 03 compute some summary tables from all the data used to generate figure panels
- one script per figure
