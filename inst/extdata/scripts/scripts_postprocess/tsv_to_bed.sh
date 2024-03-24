#!/bin/bash

# Check if input file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

input_file=$1

# Process the file
awk -F'\t' '{
    # Skip the header line
    if (NR > 1) {
        # Generate a random color for each line in the input file
        color = int(rand()*256) "," int(rand()*256) "," int(rand()*256)

        # Split the mut_max field into individual mutations
        split($4, mutations, "+")
        
        # Iterate over mutations
        for (i = 1; i <= length(mutations); i++) {
            start_field = 2*i + 3
            end_field = start_field + 1
            # Check if the start and end fields are not "NA" or process 'ref'
            if (mutations[i] == "ref" || ($start_field != "NA" && $end_field != "NA")) {
                mut_type = i":"mutations[i]  # Append index to mut_type
                # Add source element information and mutation info
                print $1"\t"$start_field"\t"$end_field"\t"mut_type"\t0\t.\t"$start_field"\t"$end_field"\t"color"\t"$2"\t"$3
            }
        }
    }
}' OFS='\t' $input_file 

