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

# Sort by reference name, then by reference start position
sort -k1,1 -k2,2n $INPUT > $TEMP

# Filtering while considering clustering of meaningless breakpoints
awk -v prox=$PROXIMITY_THRESHOLD 'BEGIN{FS=OFS="\t"}
function abs(v) {return v < 0 ? -v : v}
{
    currTName = $1; currTStart = $2; currTStrand = $7;
    currQName = $4; currQStart = $5;
    # Check if current breakpoint is close to the previous breakpoint in both reference and assembly
    if (prevTName == currTName && abs(prevTStart - currTStart) <= prox && prevQName == currQName && abs(prevQStart - currQStart) <= prox) {
        # Marking previous and current as clustered
        clustered = 1;
    } else {
        # If the previous was not part of a cluster, print it
        if (!clustered && prevLine != "") {
            print prevLine;
        }
        clustered = 0; # Reset the clustered flag
    }

    prevTName = currTName; prevTStart = currTStart; prevTStrand = currTStrand;
    prevQName = currQName; prevQStart = currQStart;

    prevLine = $0; # Store the current line as previous line for the next iteration
}
END {
    # Handle the last line
    if (!clustered) {
        print prevLine;
    }
}' $TEMP > "$OUTPUT"

# Cleanup
rm $TEMP