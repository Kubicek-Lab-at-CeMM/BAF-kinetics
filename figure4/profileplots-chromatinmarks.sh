
### profile plot chrom marks 11 clusters using deeptools

1. files=$(cat samples_chrom_features.txt)

2. for element in {cluster6,cluster7,cluster8,cluster9,cluster10,cluster11}; do matrix_1_cmd="computeMatrix reference-p oint --referencePoint center -b 1000 -a 1000 -R /data/${element}.bed -S ${files} -o /data/matrix_${element}_SL_clusters_profileplot_sorted.gz --outFileSortedRegions /data/${element}_regions_sorted_profileplot_chrom.bed" && eval $matrix_1_cmd;done

3. for element in {cluster6,cluster7,cluster8,cluster9,cluster10,cluster11}; do plotProfile -m /data/matrix_${element}_SL_clusters_profileplot_sorted.gz --perGroup -- yMax 30 --samplesLabel BRD4_Chip H3K27ac_Chip H3K4me1_Chip Pol2_Chip H3K4me3_Chip ARID1A_Chip --outFileSortedRegions /data/${element}_sorted_regions_profileplot_chrom_sorted.bed -out /data/${element}_SL_profileplot_chrom_sorted.pdf;done
