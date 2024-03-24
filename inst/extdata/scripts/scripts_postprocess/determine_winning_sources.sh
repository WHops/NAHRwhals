#!/bin/bash

# Check if merged file is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <merged_intervals_file> <winning_sources_file>"
    exit 1
fi

merged_file=$1
winning_sources_file=$2

# Determine winning sources
awk 'BEGIN { FS=OFS="\t" }
{
    # Read merged intervals and associated source information
    chr=$1; start=$2; end=$3; source_starts=$6; source_ends=$7;

    # Split the source_starts and source_ends into arrays
    split(source_starts, startArr, ",");
    split(source_ends, endArr, ",");

    # Initialize variables to track the shortest source element
    minSourceLength = "inf";
    minSourceIndex = -1;

    # Iterate through the source elements to find the shortest one that covers the entire region
    for (i = 1; i <= length(startArr); i++) {
        if ((startArr[i] <= start) && (endArr[i] >= end)) {
            sourceLength = endArr[i] - startArr[i];
            if (sourceLength < minSourceLength) {
                minSourceLength = sourceLength;
                minSourceIndex = i;
            }
        }
    }

    # Output the identifier of the winning source element and the contested region
    if (minSourceIndex > 0) {
        print chr, start, end, startArr[minSourceIndex], endArr[minSourceIndex];
    }
}' $merged_file > $winning_sources_file

echo "Winning sources determined: $winning_sources_file"
