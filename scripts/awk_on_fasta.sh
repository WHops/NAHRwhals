#!/bin/bash

# Find out the name and the length of a fasta file.

inputfa=$1
outputinfo=$2

# Run your line
awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' ${inputfa} > ${outputinfo}
