#!/bin/bash

# Semi-automatic here...

NW_res=$1


#   Res to bed including filtering and splitting rows 
./tsv_to_bed.sh <(cut -f 1,2,3,13,23,24,28,29,33,34 <(awk 'BEGIN {FS=OFS="\t"} ($11 > 98)' $NW_res))

# Sort; this was previously part of regional_dom1.sh
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$9"\t"$10"\t"$11}' output_colored_source_info.bed | sort -k1,1 -k2,2n > sorted_intervals.bed

# Preprocess for merge
./preprocess_for_merge.sh sorted_intervals.bed

# Merge
bedtools merge -i preprocessed_intervals.bed -c 4,5,6 -o distinct > merged_intervals.bed

# Postprocess
./postprocess_merged_intervals.sh merged_intervals.bed

# Winning sources
./determine_winning_sources.sh postprocessed_merged_intervals.bed

# Get the winning sources
./map_winning_sources_to_intervals.sh sorted_intervals.bed winning_sources.bed