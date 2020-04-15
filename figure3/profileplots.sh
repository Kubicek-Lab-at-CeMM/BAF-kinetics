#!/bin/bash

# deepTools needs to be installed for this to run (https://deeptools.readthedocs.io/en/develop/index.html)

# profile plot chrom marks

files=$(cat ../data/samples_chrom_features.txt)

for element in {cluster6,cluster7,cluster8,cluster9,cluster10}; do matrix_1_cmd="computeMatrix reference-point --referencePoint center -b 1000 -a 1000 -R ../data/${element}-SL-only-regions-5-cluster.txt -S ${files} -o ../data/${element}/matrix_${element}_SL_clusters_profileplot_sorted.gz --outFileSortedRegions ../data/${element}/${element}_regions_sorted_profileplot_chrom.bed" && eval $matrix_1_cmd;done

for element in {cluster6,cluster7,cluster8,cluster9,cluster10}; do plotProfile -m ../data/${element}/matrix_${element}_SL_clusters_profileplot_sorted.gz --perGroup -- yMax 25 --samplesLabel BRD4_Chip H3K27ac_Chip H3K4me1_Chip Pol2_Chip H3K4me3_Chip ARID1A_Chip --outFileSortedRegions ${element}/${element}_sorted_regions_profileplot_chrom_sorted.bed -out ${element}/${element}_SL_profileplot_chrom_sorted.pdf;done

