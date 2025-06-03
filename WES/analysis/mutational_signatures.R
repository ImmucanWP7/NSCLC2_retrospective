

# Libraries

library(reticulate)
library(SigProfilerMatrixGeneratorR)
library(signature.tools.lib)

use_virtualenv("~/Software/python/sigenv")


# Global variables

## path to folder with vcf files (one per sample)
DATA_DIR = "./data/vcf"
## output directory
OUTPUT_DIR = 'signatures_output'
## number of CPU to use for mutational signature refitting
NPARALLEL = 64



# Generate mutational catalogs (SBS)

matrices = SigProfilerMatrixGeneratorR("NSCLC2", "GRCh38", DATA_DIR, 
                                       plot=T, exome=F, bed_file=NULL, chrom_based=F, tsb_stat=F, seqInfo=F)



# Normalize mutation counts from WES to WGS
# WES kit: Twist Human Core Exome + RefSeq + Mitochondrial Panel (Twist Bioscience)

## compute normalizing factor
tricounts = read.csv("trinucleotide_counts_GRCh38.csv", stringsAsFactors=F, row.names=1)
tricounts$norm = tricounts$genome / tricounts$twist

## tansform mutation types to trinucleotides
catalog96_trinucl = sapply(rownames(matrices$'96'), function(s) paste0(substr(s,1,1),substr(s,3,3),substr(s,7,7)))

## normalize the mutational catalog by the normalizing factor
catalog96 = sweep(matrices$'96', MARGIN=1, tricounts[catalog96_trinucl,'norm'], FUN="*")

## renormalize to keep number of mutations the same
catalog96 = sweep(catalog96, MARGIN=2, colSums(matrices$'96')/colSums(catalog96), FUN=`*`)



# Refit mutational catalogs to predefined mutational signatures

## reorder rows to match COSMIC signatures
catalog96 = catalog96[match(rownames(COSMIC30_subs_signatures),rownames(catalog96)),]

## sort columns by TBS (not necessary)
catalog96 = catalog96[,order(colSums(catalog96), decreasing=T)]

## perform signature refitting using a multi-step approach where organ-specific common and rare signatures are used
sigFit = FitMS(catalogues=catalog96,
               organ='Lung',
               commonSignatureTier='T2',
               rareSignatureTier='T1',
               exposureFilterType='giniScaledThreshold',
               useBootstrap=TRUE,
               nparallel=NPARALLEL)

## save the summary plots and the results
dir.create(OUTPUT_DIR)
#plotFitMS(sigFit, outdir=OUTPUT_DIR)
#writeFitResultsToJSON(fitObj=sigFit, compress=F, filename=file.path(OUTPUT_DIR,'results_FitMS.json'))



# Save the final mutational signature decomposition to a table

## cleaning names and values
patients = sapply(strsplit(row.names(sigFit$exposures),'-'), function(x) paste0(x[c(1,2,3)],collapse='-'))
sigFit$cossim_catalogueVSreconstructed[is.na(sigFit$cossim_catalogueVSreconstructed)] = 0

## compose the final table
result = data.frame(patient=patients,
                    sample_vs_control=row.names(sigFit$exposures),
                    accuracy=sigFit$cossim_catalogueVSreconstructed,
                    sigFit$exposures,
                    row.names=NULL)

## save the resulting table
write.table(result, file=file.path(OUTPUT_DIR,'Lung_NSCLC2.tsv'), quote=F, sep='\t', row.names=F)
