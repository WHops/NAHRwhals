#!/bin/bash

# Check if input file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <sorted_intervals_file>"
    exit 1
fi

sorted_file=$1
preprocessed_file="preprocessed_intervals.bed"

# Preprocess intervals
awk 'BEGIN { FS=OFS="\t" }
{
    # Combine source start and end into a single identifier
    print $1, $2, $3, $4, $5, $6"-"$7;
}' $sorted_file > $preprocessed_file

echo "Preprocessed file created: $preprocessed_file"
