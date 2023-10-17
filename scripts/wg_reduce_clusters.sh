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

# Part A: merge clusters with 0.75 reciprocal overlap. 

# Set the number of times you want to run the command
num_runs=10

# Create a copy of the original input file
cp "$input_file" "temp_in.bed"

for ((i=1; i<=num_runs; i++)); do
    echo "Running iteration $i"

    # Run the command
    bedmap --count --echo-map-range --fraction-both 0.75 --delim '\t' <(bedtools sort -i "temp_in.bed") | sort | uniq | awk 'BEGIN {FS=OFS="\t"} {print $2,$3,$4}' > "temp_out"

    # Update the input file for the next iteration
    mv "temp_out" "temp_in.bed"
done



# Assign cluster IDs
bedtools cluster -i <(bedtools sort -i temp_in.bed) -d $cluster_dist > temp_in_2.bed

# Sort all entries by ID ($4) and subsequentlysize ($3 - $2)
awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,$3-$2,$4}' temp_in_2.bed | sort -n -r -k 5 -k 4 | cut -f 1,2,3,5 > clusters_sorted.bed

awk -v max_members=$max_members 'BEGIN {FS=OFS="\t"; prev_clusterID=-1; count=0} {if ($4 != prev_clusterID) {print $1,$2,$3,$4; prev_clusterID=$4; count=1} else if (count < max_members) {print $1,$2,$3,$4; count++}}' clusters_sorted.bed > clusters_sorted_filtered.bed

bedtools sort -i clusters_sorted_filtered.bed | cut -f 1-3 > $output_file
#
rm clusters_sorted_filtered.bed clusters_sorted.bed temp_in_2.bed temp_in.bed
