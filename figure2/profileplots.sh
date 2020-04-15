#!/bin/bash

## deepTools needs to be installed for this to run (https://deeptools.readthedocs.io/en/develop/index.html)

## profile plots tc

files=$(cat ../data/samples_tc.txt)

for element in {cluster1,cluster2,cluster3,cluster4,cluster5}; do matrix_1_cmd="computeMatrix reference-point --referencePoint center -b
1000 -a 1000 -R ../data/${element}-FC-based.bed -S ${files} -o ../data/${element}/matrix_${element}_SMARCA4_dtag_profile_plot.gz --outFileSortedRegions ../data/${element}/${element}_regions_sorted_profileplot.bed" && eval $matrix_1_cmd;done

for element in {cluster1,cluster2,cluster3,cluster4,cluster5}; do plotProfile -m ../data/${element}/matrix_${element}_SMARCA4_dtag_profile_plot.gz --yMax 25 --samplesLabel WT_SMARCA4_dtag_DMSO 2h_dtag 3h_dtag 6h_dtag 24h_dtag 72h_dtag SMARCA4_KO --perGroup -out ${element}/${element}_SMARCA4_dtag_heatmap_profile_plot.pdf;done

## profile plots chromatin features

files=$(cat samples_chrom_features.txt)

for element in {cluster1,cluster2,cluster3,cluster4,cluster5}; do matrix_1_cmd="computeMatrix reference-point --referencePoint center -b 1000 -a 1000 -R ../data/${element}-FC-based.bed -S  ${files} -o ../data/${element}/matrix_${element}_SMARCA4_dtag_profile_plot_chrom.gz --outFileSortedRegions ../data/${element}/${element}_regions_sorted_profileplot_chrom.bed" && eval $matrix_1_cmd;done

for element in {cluster1,cluster2,cluster3,cluster4,cluster5}; do plotProfile -m ../data/${element}/matrix_${element}_SMARCA4_dtag_profile_plot_chrom.gz --perGroup --yMax 25 --samplesLabel BRD4_Chip H3K27ac_Chip H3K4me1_Chip Pol2_Chip H3K4me3_Chip ARID1A_Chip --outFileSortedRegions  ${element}/${element}_sorted_regions_profileplot_chrom.bed -out ${element}/${element}_SMARCA4_dtag_profileplot_chrom.pdf;done

