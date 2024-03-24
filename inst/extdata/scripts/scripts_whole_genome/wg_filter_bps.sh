#!/bin/bash

# Check for the correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_breakpoints_file> <output_filtered_file> <proximity_threshold>"
    exit 1
fi

INPUT="$1"
OUTPUT="$2"
PROXIMITY_THRESHOLD="$3"

TEMP="temp_sorted_breakpoints.paf"
OUTPUT_PLUS="output_plus.paf"
OUTPUT_MINUS="output_minus.paf"

# Sort by reference name, then by reference start position
sort -k1,1 -k2,2n "$INPUT" > "$TEMP"

# Define the awk script as a function
process_awk() {
    awk -v prox="$PROXIMITY_THRESHOLD" 'BEGIN{FS=OFS="\t"}
    function abs(v) {return v < 0 ? -v : v}
    {
        currTName = $1; currTStart = $2; currTStrand = $7;
        currQName = $4; currQStart = $5;
        if (prevTName == currTName && abs(prevTStart - currTStart) <= prox && prevQName == currQName && abs(prevQStart - currQStart) <= prox) {
            clustered = 1;
        } else {
            if (!clustered && prevLine != "") {
                print prevLine;
            }
            clustered = 0;
        }
        prevTName = currTName; prevTStart = currTStart; prevTStrand = currTStrand;
        prevQName = currQName; prevQStart = currQStart;
        prevLine = $0;
    }
    END {
        if (!clustered) {
            print prevLine;
        }
    }' "$1"
}

# Process lines with "+" in $7
process_awk <(grep -e '\+$' "$TEMP") > "$OUTPUT.plus"

# Process lines with "-" in $7
process_awk <(grep -e '-$' "$TEMP") > "$OUTPUT.minus"

# Concatenate the results and sort them
cat "$OUTPUT.plus" "$OUTPUT.minus" | sort -k1,1 -k2,2n > "$OUTPUT"

# Cleanup
rm "$TEMP" "$OUTPUT.plus" "$OUTPUT.minus" # Remove temporary files if not needed