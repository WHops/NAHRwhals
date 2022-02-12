#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified March 8, 2017

Description:  Concatenates sam files, keeping only the header from the first.

Usage:  mergesam.sh <files> out=<file>

Java Parameters:
-da             Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

#This block allows symlinked shellscripts to correctly set classpath.
pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx400m"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	setEnvironment
	parseXmx "$@"
}
calcXmx "$@"

function mergesam() {
	local CMD="java $EA $EOOM $z -cp $CP jgi.MergeSam $@"
	echo $CMD >&2
	eval $CMD
}

mergesam "$@"
