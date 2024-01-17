#!/bin/bash

# Check for the correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: ./bedtools_per_querycontig.sh INPUT OUTPUT MERGE_DISTANCE"
    exit 1
fi

INPUT="$1"
OUTPUT="$2"
MERGE_DISTANCE="$3"
TEMP="temp_reordered.bed"
TEMP2="temp_merged.bed"

# Rearrange columns and sort
awk 'OFS="\t" {print $4"_"$1, $2, $3, $1}' "$INPUT" | sort -k1,1 -k2,2n > "$TEMP"

# Merge using bedtools
bedtools merge -i "$TEMP" -d $MERGE_DISTANCE -c 4,1 -o distinct > "$TEMP2"

# Restore original order and format
awk '{FS=OFS="\t"} ($3-$2 > 2) {print $4, $2, $3, $1}' "$TEMP2" > "$OUTPUT"

# Cleanup
rm "$TEMP" "$TEMP2"