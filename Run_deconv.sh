#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd $SCRIPT_DIR
ref_path=$1
val_path=$2
props_path=$3
p_val_cutoff=$4
n_markers=$5

Rscript ./Scripts/IC_deconvolution_workflow_part1.r $ref_path $val_path $props_path $p_val_cutoff $n_markers
python ./Scripts/PerformMethatlas.py $props_path
Rscript ./Scripts/IC_deconvolution_workflow_part2.r $props_path