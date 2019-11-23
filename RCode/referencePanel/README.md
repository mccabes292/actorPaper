# R code for using GTEx as a Reference Panel

* `getColNames.R` - Calculates the column names of the isoform matrix for easy subsetting of the demographic files
* `calcLibrarySize.R` - Calculates library size of GTEx
* `calcScaleTPM.R` - Creates scaled TPM matrix for GTEx 
* `calcPrecisionTPM_NoIsoFilt_MLE.R` - Takes as input an index corresponding a GTEx tissue and precomputes Dirichlet Multinomial maximum likelihood estimates for each gene. Writes out results per tissue.
* `combinePrecisionNoIsoFilt_MLE.R` - Combines results across all tissues from `calcPrecisionTPM_NoIsoFilt_MLE.R`
* `reduceAlphaDir.R` - Reduces the combined Dirichlet Multinomial estimates to only genes with less than 6 annotated isoforms
