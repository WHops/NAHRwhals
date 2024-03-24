#!/bin/bash

# Semi-automatic here...

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <NW_res> <OUTFILE> <OUTDIR> <SCRIPTDIR_BASE>"
    exit 1
fi

NW_res=$1
OUTFILE=$2
OUTDIR=$3
SCRIPTDIR_BASE=$4

mkdir -p $OUTDIR

# In the input file, I want to check which ones have $11 < 98. For those that do, then column 13 should become 'ref'
awk 'BEGIN {FS=OFS="\t"} ($11 < 98) {$13="ref"; print $0}' $NW_res > $OUTDIR/res.ref

# Next, wherever $13 is ref, now I want $23 and $24 to be the same as col 2 and 3, respectively 
awk 'BEGIN {FS=OFS="\t"} ($13 == "ref") {$23=$2; $24=$3 ; print $0}' $OUTDIR/res.ref > $OUTDIR/res.ref2


#  Res to bed including filtering and splitting rows 
$SCRIPTDIR_BASE/tsv_to_bed.sh <(cut -f 1,2,3,13,23,24,28,29,33,34 <(awk 'BEGIN {FS=OFS="\t"} (($11 > 98) && ($13 != "ref"))' $NW_res)) > ${OUTDIR}/output_colored_source_info.bed

# make here an exit: if ${OUTDIR}/output_colored_source_info.bed is empty, print a warning and exit.
if [ ! -s ${OUTDIR}/output_colored_source_info.bed ]; then
    echo "No SVs detected! Exiting."
    exit 1
fi

# Sort; this was previously part of regional_dom1.sh
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$9"\t"$10"\t"$11}' ${OUTDIR}/output_colored_source_info.bed | sort -k1,1 -k2,2n > ${OUTDIR}/sorted_intervals.bed

# Preprocess for merge
$SCRIPTDIR_BASE/preprocess_for_merge.sh ${OUTDIR}/sorted_intervals.bed ${OUTDIR}/preprocessed_intervals.bed

# Merge
bedtools merge -i ${OUTDIR}/preprocessed_intervals.bed -c 4,5,6 -o distinct > ${OUTDIR}/merged_intervals.bed

# Postprocess
$SCRIPTDIR_BASE/postprocess_merged_intervals.sh ${OUTDIR}/merged_intervals.bed ${OUTDIR}/postprocessed_merged_intervals.bed

# Add ref associations to the merged invervals
$SCRIPTDIR_BASE/memory_inception_refsources_to_merges.sh ${OUTDIR}/postprocessed_merged_intervals.bed $OUTDIR/res.ref2 $OUTDIR/updated_postprocessed.bed

# Winning sources
$SCRIPTDIR_BASE/determine_winning_sources.sh $OUTDIR/updated_postprocessed.bed $OUTDIR/winning_sources.bed

# Get the winning sources
$SCRIPTDIR_BASE/map_winning_sources_to_intervals.sh $OUTDIR/sorted_intervals.bed $OUTDIR/winning_sources.bed $OUTDIR/dominant_intervals.bed

$SCRIPTDIR_BASE/merge_final_stuff.sh $OUTDIR/dominant_intervals.bed $OUTFILE


