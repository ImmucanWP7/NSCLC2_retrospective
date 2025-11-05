# NSCLC2_retrospective

This repo contains the code for the publication of the analysis of the retrospective NSCLC2 cohort from IMMUcan.
Most folders are set up in a `workflowR` manner.

# Data

The data is available upon request from the IMMUcan consortium. Further details will be provided

# Folders

## IMC

1. read IMC panel 1 data
2. perform cell typing

## IMC_p2

1. read IMC panel 1 data
2. perform cell typing

## RNA

1. prepare metadata
2. RNA count processing 
3. QC plots
4. deconvolution or RNAseq data

## WES

- immucan_nsclc_wes_tableReport.Rmd mutation identification and Oncoplot
- mutational_signatures.R generation of mutational signatures

## MOFA

1. data preparation
2. model training

## Figures

contains scripts to generate Figures. Depends on all other folders to be run first.
