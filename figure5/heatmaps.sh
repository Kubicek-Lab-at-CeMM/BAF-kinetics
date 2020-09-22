#!/bin/bash

## deepTools needs to be installed for this to run (https://deeptools.readthedocs.io/en/develop/index.html)

files=$(cat ../data/samples_late_timepoints_withARID2KO.txt)

matrix_1_cmd="computeMatrix reference-point --referencePoint center -b 1500 -a 1500 -R ../data/supenh-diff-any-comp.bed enhancer-diff-wo-supenhancer-all.bed -S ${files} -o ../data/matrix.H3K27ac.broad.wlatetimepoints.gz --outFileSortedRegions ../data/regions_sorted_H3K27ac_broad_wlatetimepoints.bed"

eval $matrix_1_cmd

plotHeatmap -m ../data/matrix.H3K27ac.broad.wlatetimepoints.gz --sortUsingSamples 10 --zMax 25 --yMax 35 --missingDataColor 1 --colorMap 'Blues' --samplesLabel WT_SMARCA4_dtag_DMSO 72h_dtag A2_A4_F5_DMSO A2_A4_F5_dTAG_47_72h BRG_BRM_DMSO 72h_dtag CC1_CC2_DMSO 72h_dtag BRD4_Chip H3K27ac_Chip H3K4me1_Chip ARID1A_Chip --outFileSortedRegions regions_sortedH3K27ac_after_plot_broad_H3K27acsorted.justlatetimepoints.txt -out heatmap_enh_supenh_sorted_H3K27ac_0white_1500bp_justlatetimepoints.pdf
