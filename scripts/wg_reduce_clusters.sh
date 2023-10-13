#!/bin/bash

# Check if the correct number of arguments is provided
if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <input_file> <output_file> <max_cluster_members> <max_distance_to_same_cluster"
    exit 1
fi

input_file=$1
output_file=$2
max_members=$3
cluster_dist=$4
# Assign cluster IDs
bedtools cluster -i $input_file -d $cluster_dist > clustered.bed

# Sort all entries by ID ($4) and subsequentlysize ($3 - $2)
awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,$3-$2,$4}' clustered.bed | sort -n -r -k 5 -k 4 | cut -f 1,2,3,5 > clusters_sorted.bed

awk -v max_members=$max_members 'BEGIN {FS=OFS="\t"; prev_clusterID=-1; count=0} {if ($4 != prev_clusterID) {print $1,$2,$3,$4; prev_clusterID=$4; count=1} else if (count < max_members) {print $1,$2,$3,$4; count++}}' clusters_sorted.bed > final_final.bed

bedtools sort -i final_final.bed | cut -f 1-3 > output_preoverlapcut.bed $output_file

bedmap --count --echo-map-range --fraction-both 0.9 --delim '\t' output_preoverlapcut.bed \
    | awk '$1>1' - \
    | cut -f2- - \
    | sort-bed - \
    | uniq - \
    > mutuallyOverlappingIntervals.bed

bedops --merge output_preoverlapcut.bed > mergedIntervals.bed
bedmap --echo-map --exact --skip-unmapped output_preoverlapcut.bed mergedIntervals.bed > exclusiveIntervals.bed
bedops --everything exclusiveIntervals.bed mutuallyOverlappingIntervals.bed > $output_file


rm final_final.bed clusters_sorted.bed clustered.bed

echo "Processing completed. Results saved in $output_file."