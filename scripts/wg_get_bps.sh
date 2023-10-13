#!/bin/bash

# Check for the correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_paf_file> <output_file> <skip_starts_bool>"
    exit 1
fi

INPUT="$1"
OUTPUT="$2"
SKIPSTARTS="$3"
TEMP="temp_filtered_length.paf"
TEMP2="temp_count_contigs.paf"
TEMP3="temp_breakpoints.paf"
echo "Input file: $INPUT"
# Distance from the start or end of contig to consider a breakpoint as edge
EDGE_DISTANCE=50

# Step 1: Filter out alignments shorter than 10 kbp based on the query sequence length
echo "Step 1: Filtering alignments shorter than 10 kbp based on the query sequence length..."
awk 'BEGIN{FS=OFS="\t"} ($4-$3 >= 1000) {print $0}' "$INPUT" > $TEMP
echo "Filtering done."

# Step 2: Count occurrences of each contig
echo "Step 2: Counting occurrences of each contig..."
awk 'BEGIN{FS=OFS="\t"} {contig_count[$1]++} END {for (contig in contig_count) print contig, contig_count[contig]}' $TEMP > $TEMP2
echo "Counting done."

# Step 3: Filter contigs that appear in at least 3 separate alignments and extract breakpoints
echo "Step 3: Filtering contigs with less than 3 alignments and extracting breakpoints..."
awk -v edge="$EDGE_DISTANCE" -v skipstarts="$SKIPSTARTS" 'BEGIN{FS=OFS="\t"} 
NR==FNR {contigs[$1]=$2; next} contigs[$1] >= 3 {
    qName = $1; qStart = $3; qEnd = $4; qLen = $2; strand = $5;
    tName = $6; tStart = $8; tEnd = $9;

    # Swap qstart and qend if on the negative strand
    if (strand == "-") {
        temp = qStart;
        qStart = qEnd;
        qEnd = temp;
    }

    # Skip breakpoints at the very start or end of a contig based on the skipstarts flag
    if (skipstarts) {
        if (qStart > edge && qStart < qLen - edge) {
            print tName, tStart, tStart+1, qName, qStart, qStart+1, strand;
        }
        if (qEnd > edge && qEnd < qLen - edge) {
            print tName, tEnd, tEnd+1, qName, qEnd, qEnd+1, strand;
        }
    } else {
        print tName, tStart, tStart+1, qName, qStart, qStart+1, strand;
        print tName, tEnd, tEnd+1, qName, qEnd, qEnd+1, strand;
    }
}' "$TEMP2" "$TEMP" > "$TEMP3"
echo "Extraction done."

echo "Moving the results to the desired output file..."
mv $TEMP3 "$OUTPUT"

# Cleanup
rm $TEMP $TEMP2

echo "All done! Results saved to $OUTPUT."