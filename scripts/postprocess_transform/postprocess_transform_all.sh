#!/bin/bash

# Semi-automatic here...

NW_res=$1


# In the input file, I want to check which ones have $11 < 98. For those that do, then column 13 should become 'ref'
awk 'BEGIN {FS=OFS="\t"} ($11 < 98) {$13="ref"; print $0}' $NW_res > $NW_res.ref

# Next, wherever $13 is ref, now I want $23 and $24 to be the same as col 2 and 3, respectively 
awk 'BEGIN {FS=OFS="\t"} ($13 == "ref") {$23=$2; $24=$3 ; print $0}' $NW_res.ref > $NW_res.ref2


#  Res to bed including filtering and splitting rows 
# ./tsv_to_bed.sh <(cut -f 1,2,3,13,23,24,28,29,33,34 $NW_res.ref2)
./tsv_to_bed.sh <(cut -f 1,2,3,13,23,24,28,29,33,34 <(awk 'BEGIN {FS=OFS="\t"} (($11 > 98) && ($13 != "ref"))' $NW_res))

# Sort; this was previously part of regional_dom1.sh
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$9"\t"$10"\t"$11}' output_colored_source_info.bed | sort -k1,1 -k2,2n > sorted_intervals.bed

# Preprocess for merge
./preprocess_for_merge.sh sorted_intervals.bed

# Merge
bedtools merge -i preprocessed_intervals.bed -c 4,5,6 -o distinct > merged_intervals.bed

# Postprocess
./postprocess_merged_intervals.sh merged_intervals.bed

# Add ref associations to the merged invervals
./memory_inception_refsources_to_merges.sh postprocessed_merged_intervals.bed $NW_res.ref2

# Winning sources
./determine_winning_sources.sh updated_postprocessed.bed

# Get the winning sources
./map_winning_sources_to_intervals.sh sorted_intervals.bed winning_sources.bed