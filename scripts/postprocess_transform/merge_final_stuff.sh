#!/bin/bash

input_file=$1 
output_file="merged_intervals.bed"

# Sort the file by chromosome and start position
sort -k1,1 -k2,2n "$input_file" > sorted_input.bed

# Use awk to process the sorted file
awk 'BEGIN {OFS="\t"}
    function flush_interval() {
        if (interval_start > 0 && interval_end > 0) {
            name = substr(names, 2) # Remove leading '+'
            print chr, interval_start, interval_end, name
        }
    }

    function add_name(name) {
        gsub(/[0-9:_]/, "", name) # Remove numbers, colons, and underscores
        if (!(name in name_set)) {
            names = names "+" name
            name_set[name]
        }
    }

    {
        current_chr = $1
        current_start = $2
        current_end = $3
        current_name = $4

        if (current_chr == chr && current_start <= interval_end && current_end > interval_start && !(current_start == interval_end || current_end == interval_start)) {
            # Merge with current interval
            interval_end = (current_end > interval_end ? current_end : interval_end)
            add_name(current_name)
        } else {
            # Flush current interval and start a new one
            flush_interval()
            chr = current_chr
            interval_start = current_start
            interval_end = current_end
            delete name_set
            names = ""
            add_name(current_name)
        }
    }
    END {
        flush_interval()
    }
' sorted_input.bed > "$output_file"

# Cleanup
rm sorted_input.bed

echo "Merged intervals have been saved to $output_file"