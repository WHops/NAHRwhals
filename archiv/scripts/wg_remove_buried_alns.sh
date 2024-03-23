#!/bin/bash

# Check for the correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_paf_file> <output_file>"
    exit 1
fi

INPUT="$1"
OUTPUT="$2"
TEMP_SORTED="temp_sorted.paf"

# Sort the PAF file by tName, tStart and descending tEnd
sort -k6,6 -k8,8n -k9,9nr $INPUT > $TEMP_SORTED

# Filtering contained alignments
awk 'BEGIN{FS=OFS="\t"}
{
    if ($6 != prevTName || $8 < prevTStart || $9 > prevTEnd) {
        print $0;
        prevTName = $6;
        prevTStart = $8;
        prevTEnd = $9;
    }
    # Otherwise, the alignment is contained within the previous and will be skipped
}' $TEMP_SORTED > $OUTPUT

# Cleanup
rm $TEMP_SORTED

echo "Filtered PAF saved to $OUTPUT."
