#!/bin/bash

# Check if merged file is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <merged_intervals_file> <postprocessed_file>"
    exit 1
fi

merged_file=$1
postprocessed_file=$2

# Process the merged intervals file
awk 'BEGIN { FS=OFS="\t" }
{
    chr=$1; start=$2; end=$3; mutations=$4; colors=$5; combined_sources=$6;
    split(combined_sources, sources, ",");

    source_starts = source_ends = "";

    # Process each combined source and split them
    for (i = 1; i <= length(sources); i++) {
        split(sources[i], src, "-");
        source_starts = source_starts src[1] (i < length(sources) ? "," : "");
        source_ends = source_ends src[2] (i < length(sources) ? "," : "");
    }

    # Print the postprocessed line
    print chr, start, end, mutations, colors, source_starts, source_ends;
}' $merged_file > $postprocessed_file

echo "Postprocessed intervals file created: $postprocessed_file"
