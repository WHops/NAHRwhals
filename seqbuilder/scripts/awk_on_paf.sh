#!/bin/bash

# Run awk to correct paf if it's made on scattered input

inputpaf=$1
outputpaf=$2

# Run your line
awk '{FS=OFS="\t"} match($1, /_([0-9]*)?-/,a){$3 = $3+a[1]; $4 = $4+a[1]; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' ${inputpaf} > ${outputpaf}
