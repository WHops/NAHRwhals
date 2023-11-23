#!/bin/bash

# Check if input file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

input_file=$1
output_file="output_colored_source_info.bed"

# Process the file
awk -F'\t' '{
    # Skip the header line and lines with "ref" in the mut_max column
    if (NR > 1 && $4 != "ref") {
        # Generate a random color for each line in the input file
        color = int(rand()*256) "," int(rand()*256) "," int(rand()*256)

        # Iterate over mutations
        for (i = 1; i <= 3; i++) {
            start_field = 2*i + 3
            end_field = start_field + 1
            # Check if the start and end fields are not "NA"
            if ($start_field != "NA" && $end_field != "NA") {
                mut_type = $4  # mut_max column
                # Add source element information and mutation info
                print $1"\t"$start_field"\t"$end_field"\t"mut_type"\t0\t.\t"$start_field"\t"$end_field"\t"color"\t"$2"\t"$3
            }
        }
    }
}' OFS='\t' $input_file > $output_file

echo "BED file created: $output_file"
