#!/bin/bash

# Visualize two fasta sequences fast. 

queryfasta=$1
targetfasta=$2
outname=$3
del=$4

dotplotlyscript="/Users/hoeps/PhD/projects/nahrcall/nahrchainer/seqbuilder/R/pafCoordsDotPlotly.R"

# Align
minimap2 -x asm20 \
         -c \
         -z400,50 \
         -s 0 \
         -M 0.2 \
         -N 100 \
         -P \
          ${queryfasta} \
          ${targetfasta}  > ${outname}.paf

# Plot
${dotplotlyscript} -i ${outname}.paf \
                   -o ${outname} \
                   -m 50 -p 10 -q 10

# Open
open ${outname}.pdf

# Optionally, remove traces
if [ $del -ge 1 ]; then
    rm -r ${outname}.paf ${outname}.pdf ${outname}.html ${outname}_files 
fi
