#!/bin/bash

# Check if sorted and winning sources files are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <sorted_intervals_file> <winning_sources_file>"
    exit 1
fi

sorted_file=$1
winning_sources_file=$2
output_file="dominant_intervals.bed"

# Process the files
awk 'BEGIN {OFS="\t"} NR==FNR {
    # Store winning source information with their corresponding merged intervals
    # Format: source_start_source_end -> merged_start1_merged_end1 merged_start2_merged_end2 ...
    winners[$4"_"$5] = winners[$4"_"$5] " " $2"_"$3;
    next;
}
{
    key = $6"_"$7; # Source start and end from sorted_intervals
    if (key in winners) {
        n = split(winners[key], merged_intervals, " ");
        for (i = 1; i <= n; i++) {
            split(merged_intervals[i], merged, "_");
            merged_start = merged[1];
            merged_end = merged[2];
            if ($2 >= merged_start && $3 <= merged_end) {
                # If the sorted interval is encapsulated by any of the merged regions
                print $1, $2, $3, $4, $5, ".", $2, $3, $9, $6, $7;
                break;
            }
        }
    }
}' $winning_sources_file $sorted_file > $output_file

echo "Final dominant intervals processed: $output_file"
