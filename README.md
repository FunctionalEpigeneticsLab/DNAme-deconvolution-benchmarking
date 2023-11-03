# DNAme-deconvolution-benchmarking

The scripts provide details for the deconvolutions performed in the manuscript entitled 'Benchmarking of Methods for DNA Methylome Deconvolution'. These include Python, Bash and R scripts.

> This is not the release of a software package. All scripts and binaries are provided as is, without any warranty. We are only providing this information and code in addition to the description of methods for making it easier to reproduce deconvolutions.

## Versions
```
Python 3.8.11
R version 4.1.1
```
## Dependencies
```
R-packages:
  - glmnet
  - bvls
  - preprocessCore
  - quadprog
  - stringr
  - nnls
  - EpiDISH
  - minfi
  - EMeth
  - Matrix
  - MethylResolver
  - TOAST
  - corpcor
  - writexl
  - limma

Python-packages:
  - numpy
  - pandas
  - scipy
  - multiprocessing
```


## Bioinformatics processing

Sequencing reads were trimmed with TrimGalore! (v0.6.6) to remove potential adapter contamination and to trim off the adaptase-induced addition of random nucleotides. 

Next, trimmed reads were aligned to the reference human genome GRCh37 using Bwameth (v0.2.5). Reads were deduplicated using Picard (v3.1.0) and methylation was called using MethylDackel (v0.6.0). The pipeline used for this can be found in the following repo: https://github.com/KobeDR/Bwameth.

Further analyses were performed using R (v4.1.1) and Python (v3.8.11.).

## Instructions
Input arguments:
  - Path to reference data frame: csv format file containing methylation ratios of reference dataset; rows depicting loci, columns depicting biological replicates of cell types (names of cell types should occur in the column names e.g., 'Neu.1' is valid if 'Neu' is a cell type in the proportions file).
  - Path to validation data frame: csv format file containing methylation ratios of validation dataset; rows depicting loci, columns depicting mixtures.
  - Path to respective cell proportions in validation data frame:  csv format file containing proportions of validation dataset; rows depicting mixtures, columns depicting cell types.
  - P value cutoff for marker selection (recommended: 0.05).
  - Preferred number of markers per cell type (recommended: 100).

Output argument:
  - Csv file containing predicted proportions along with actual proportions, predicted cell type and normalization and deconvolution algorithm that was applied.

Example code:
```
  ./Run_deconv.sh /path/to/ref_df.csv /path/to/val_df.csv /path/to/props_df.csv 0.05 100
```
