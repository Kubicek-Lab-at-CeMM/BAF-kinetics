#!/bin/bash

## deepTools needs to be installed for this to run (https://deeptools.readthedocs.io/en/develop/index.html)

files=$(cat ../data/samples_woARID2KO.txt)

matrix_1_cmd="computeMatrix reference-point --referencePoint center -b 1000 -a 1000 -R ../data/differential-enhancer-bedfiles.txt -S ${files} -o ../data/Enhancer_matrix.gz --outFileSortedRegions ../data/Enhancer_regions_sorted.bed"

eval($matrix_1_cmd)

plotHeatmap -m ../data/Enhancer_matrix.gz --zMax 25 --yMax 25 --sortUsingSamples 24 --colorMap 'Blues' --missingDataColor 1 --samplesLabel WT_SMARCA4_dtag_DMSO 2h_dtag 3h_dtag 6h_dtag 24h_dtag 72h_dtag SMARCA4_KO SMARCA2_KO_DMSO BRG_BRM_DMSO 2h_dtag 3h_dtag 6h_dtag  24h_dtag 72h_dtag SMARCC2_KO_DMSO SMARCC1_KO_DMSO CC1_CC2_DMSO 2h_dtag 3h_dtag 6h_dtag 24h_dtag 72h_dtag BRD4_Chip H3K27ac_Chip H3K4me1_Chip Pol2_Chip H3K4me3_Chip ARID1A_Chip SMARCA4_Chip --outFileSortedRegions Enhancer_regions_sorted_plot.bed -out Enhancer_Heatmap_sorted_0white.pdf
