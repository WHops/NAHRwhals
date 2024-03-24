#!/bin/bash

# A script to process genomic data. 
# This script includes various tools for data manipulation such as minimap2, bedtools, and custom scripts.

# Use: ./script_name.sh [input_file] [output_dir] [genome_path] [bedtools_path]
# Example: ./script_name.sh wga.paf ~/results/ ~/PhD/projects/huminvs/genomes/hg38/hg38.genome /path/to/bedtools

set -euo pipefail

# Check if correct number of arguments are provided
if [ "$#" -ne 9 ]; then
    echo "Usage: $0 [input_file] [output_directory] [output_file] [genome_path] [bedtools_path] [cluster_merge_distance] [indel_ignore_threshold] [exclusion_mask] [scripts_basedir]"
    exit 1
fi

INPUT="$1"
OUTPUT_DIR="$2"
OUTPUT_FILE="$3"
GENOME_PATH="$4"
BEDTOOLS="$5"
BP_MERGE_DISTANCE="$6"
BP_COLLAPSE_REMOVE_DISTANCE="$7"
EXCLUSION_MASK="$8"
SCRIPTDIR_BASE="$9"
EXPAND_FRACTION=1
MERGE_DISTANCE=1000000
FINAL_REGIONS_CLUSTER_DISTANCE=100000
FINAL_INTERVALS_PER_CLUSTER=3

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Step 1: Remove buried alignments
${SCRIPTDIR_BASE}/wg_remove_buried_alns.sh "$INPUT" "${OUTPUT_DIR}/wga_primary_nonoverlap.paf"

# Step 2: Extract breakpoints
${SCRIPTDIR_BASE}/wg_get_bps.sh "${OUTPUT_DIR}/wga_primary_nonoverlap.paf" "${OUTPUT_DIR}/wga_primary_bps.bed" 1
${SCRIPTDIR_BASE}/wg_get_bps.sh "${OUTPUT_DIR}/wga_primary_nonoverlap.paf" "${OUTPUT_DIR}/wga_primary_bps_noskip.bed" 0

# Step 3: Filter breakpoints based on a distance of 10,000
${SCRIPTDIR_BASE}/wg_filter_bps.sh "${OUTPUT_DIR}/wga_primary_bps.bed" "${OUTPUT_DIR}/wga_primary_bps_distfiltered.bed" $BP_COLLAPSE_REMOVE_DISTANCE

# Step 4: Extract clusters based on a distance of 1,000,000
${SCRIPTDIR_BASE}/wg_extract_clusters.sh "${OUTPUT_DIR}/wga_primary_bps_distfiltered.bed" "${OUTPUT_DIR}/breakpoint_clusters_merged.bed" $BP_MERGE_DISTANCE

# Step 5: Expand the clusters
$BEDTOOLS slop -i "${OUTPUT_DIR}/breakpoint_clusters_merged.bed" -g "$GENOME_PATH" -b $EXPAND_FRACTION -pct > "${OUTPUT_DIR}/breakpoint_clusters_merged_expand.bed"

# Step 6: Filter clusters greater than 50,000 after expansion
cp "${OUTPUT_DIR}/breakpoint_clusters_merged_expand.bed" "${OUTPUT_DIR}/breakpoint_clusters_merged_expand_50kb.bed"

# Step 7: Filter out random, alt, and chrUn entries, and then sort
grep -vE 'random|alt|chrUn' "${OUTPUT_DIR}/breakpoint_clusters_merged_expand_50kb.bed" | $BEDTOOLS sort -i - | awk 'BEGIN {FS=OFS="\t"} {print $4,$2,$3,$1}' > "${OUTPUT_DIR}/list.bed"
grep -vE 'random|alt|chrUn' "${OUTPUT_DIR}/breakpoint_clusters_merged_expand_50kb.bed" | $BEDTOOLS sort -i - | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3}' > "${OUTPUT_DIR}/list_igv.bed"


# Step 8: Formatting, sorting, merging, and again formatting the data
awk 'BEGIN {FS=OFS="\t"} {print $1"_"$6,$8,$9,$6}' "${OUTPUT_DIR}/wga_primary_nonoverlap.paf" | $BEDTOOLS sort -i - | $BEDTOOLS merge -d $MERGE_DISTANCE -c 4 -o distinct -i - | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,$4}' > "${OUTPUT_DIR}/contigmerges.bed"
awk 'BEGIN {FS=OFS="\t"} {print $1"_"$6,$8,$9,$6}' "${OUTPUT_DIR}/wga_primary_nonoverlap.paf" | $BEDTOOLS sort -i - | $BEDTOOLS merge -d $MERGE_DISTANCE -c 4 -o distinct -i - | awk 'BEGIN {FS=OFS="\t"} {print $4,$2,$3,$1}' > "${OUTPUT_DIR}/contigmerges_igv.bed"

# Step 9: Intersect the processed data
$BEDTOOLS intersect -a "${OUTPUT_DIR}/contigmerges.bed" -b "${OUTPUT_DIR}/list.bed" | awk 'BEGIN {FS=OFS="\t"} {print $4,$2,$3}' > "${OUTPUT_DIR}/list_cut_final_prefilter.bed"

${SCRIPTDIR_BASE}/wg_reduce_clusters.sh <(bedtools sort -i ${OUTPUT_DIR}/list_cut_final_prefilter.bed) ${OUTPUT_DIR}/list_cut_final.bed $FINAL_INTERVALS_PER_CLUSTER $FINAL_REGIONS_CLUSTER_DISTANCE

if [ "$EXCLUSION_MASK" != "none" ]; then
    $BEDTOOLS subtract -a "${OUTPUT_DIR}/list_cut_final.bed" -b "$EXCLUSION_MASK" > "${OUTPUT_DIR}/list_cut_final_exclusion.bed"
    mv "${OUTPUT_DIR}/list_cut_final_exclusion.bed" $OUTPUT_FILE
fi

