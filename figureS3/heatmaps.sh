#!/bin/bash

## deepTools needs to be installed for this to run (https://deeptools.readthedocs.io/en/develop/index.html)

files=$(cat ../data/samples.txt)

for element in {cluster1,cluster2,cluster3,cluster4,cluster5}; do matrix_1_cmd="computeMatrix reference-point --referencePoint center -b 1000 -a 1000 -R ../data/${element}-FC-based.bed -S ${files} -o ../data/${element}/matrix_${element}_SMARCA4_dtag.gz --outFileSortedRegions ../data/${element}/${element}_regions_sorted.bed" && eval $matrix_1_cmd;done

for element in {cluster1,cluster2,cluster3,cluster4,cluster5}; do plotHeatmap -m ../data/${element}/matrix_${element}_SMARCA4_dtag.gz --zMax 25 --yMax 25 --samplesLabel WT_DMSO WT_SMARCA4_dtag_DMSO 2h_dtag 3h_dtag 6h_dtag 24h_dtag 72h_dtag SMARCA4_KO BRD4_Chip H3K27ac_Chip H3K4me1_Chip Pol2_Chip H3K4me3_Chip ARID1A_Chip SMARCA4_Chip --outFileSortedRegions ${element}/${element}_sorted_regions.bed -out ${element}/${element}_SMARCA4_dtag_heatmap.png;done
