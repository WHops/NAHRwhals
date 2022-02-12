#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified July 10, 2017

Description:  Tests the memory usage of a sequence file.

Usage:  loadreads.sh in=<file>

Parameters:
in=             Input file.
lowcomplexity=f Assume input library is low-complexity.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify
                20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max
                is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
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

z="-Xmx2g"
z2="-Xms2g"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	setEnvironment
	parseXmx "$@"	
	if [[ $set == 1 ]]; then
		return
	fi
	freeRam 2000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

testmem() {
	local CMD="java $EA $EOOM $z $z2 -cp $CP driver.LoadReads $@"
	echo $CMD >&2
	eval $CMD
}

testmem "$@"
