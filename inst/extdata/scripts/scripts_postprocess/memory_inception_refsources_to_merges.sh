#!/bin/bash

# Define the file names
postprocessed_file=$1
sorted_file=$2
output_file=$3

# Create or clear the output file
> "$output_file"

# Read each line from the postprocessed file
while IFS=$'\t' read -r chr start end name color assoc_starts assoc_ends; do
    # Initialize variables to store updated associated starts and ends
    new_assoc_starts="$assoc_starts"
    new_assoc_ends="$assoc_ends"

    # Split the associated starts and ends into arrays
    IFS=',' read -ra starts_array <<< "$assoc_starts"
    IFS=',' read -ra ends_array <<< "$assoc_ends"

    # Read each line from the sorted file
    while IFS=$'\t' read -r sorted_chr sorted_start sorted_end sorted_name sorted_irrelevant1 sorted_irrelevant2; do
        # Check for overlap and the same chromosome
        if [[ "$chr" == "$sorted_chr" ]] && ((sorted_start < end)) && ((sorted_end > start)); then
            # Check if sorted interval is fully contained within any associated intervals
            for i in "${!starts_array[@]}"; do
                if ((sorted_start >= starts_array[i])) && ((sorted_end <= ends_array[i])); then
                    # Append to the associated starts and ends
                    new_assoc_starts="$new_assoc_starts,$sorted_start"
                    new_assoc_ends="$new_assoc_ends,$sorted_end"
                    break # Break the loop as we found a containing interval
                fi
            done
        fi
    done < "$sorted_file"

    # Write the updated line to the output file
    echo -e "$chr\t$start\t$end\t$name\t$color\t$new_assoc_starts\t$new_assoc_ends" >> "$output_file"
done < "$postprocessed_file"

echo "Updated file is saved as $output_file"