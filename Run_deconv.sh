#!/bin/bash
#### Important! If software dependencies have not been installed yet, run './Install_depencies.sh'
#### Pass the following arguments in this order:
#   - Path to data frame of reference methylation values (csv format - rows depicting CpGs - columns depicting samples).
#   - Path to data frame of validation methylation values (csv format - rows depicting CpGs - columns depicting samples).
#   - Path to data frame of actual proportions (csv format - rows depicting samples - columns depicting cell types).
ref_path=$1
val_path=$2
props_path=$3
Rscript ./Scripts/IC_deconvolution_workflow_part1.r $ref_path $val_path
python ./Scripts/PerformMethatlas.py $props_path
Rscript ./Scripts/IC_deconvolution_workflow_part2.r $props_path