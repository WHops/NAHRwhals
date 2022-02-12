#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified Jan 7, 2020

Description:  Demultiplexes sequences into multiple files based on their names,
substrings of their names, or prefixes or suffixes of their names.
Allows unlimited output files while maintaining only a small number of open file handles.

Usage:
demuxbyname.sh in=<file> in2=<file2> out=<outfile> out2=<outfile2> names=<string,string,string...>

Alternately:
demuxbyname.sh in=<file> out=<outfile> delimiter=whitespace prefixmode=f
This will demultiplex by the substring after the last whitespace.

demuxbyname.sh in=<file> out=<outfile> length=8 prefixmode=t
This will demultiplex by the first 8 characters of read names.

demuxbyname.sh in=<file> out=<outfile> delimiter=: prefixmode=f
This will split on colons, and use the last substring as the name; useful for
demuxing by barcode for Illumina headers in this format:
@A00178:73:HH7H3DSXX:4:1101:13666:1047 1:N:0:ACGTTGGT+TGACGCAT

in2 and out2 are for paired reads in twin files and are optional.
If input is paired and there is only one output file, it will be written interleaved.

File Parameters:
in=<file>       Input file.
in2=<file>      If input reads are paired in twin files, use in2 for the second file.
out=<file>      Output files for reads with matched headers (must contain % symbol).
                For example, out=out_%.fq with names XX and YY would create out_XX.fq and out_YY.fq.
                If twin files for paired reads are desired, use the # symbol.  For example,
                out=out_%_#.fq in this case would create out_XX_1.fq, out_XX_2.fq, out_YY_1.fq, etc.
outu=<file>     Output file for reads with unmatched headers.
stats=<file>    Print statistics about how many reads went to each file.

Processing Modes (determines how to convert a read into a name):
prefixmode=t    (pm) Match prefix of read header.  If false, match suffix of read header.
                prefixmode=f is equivalent to suffixmode=t.
barcode=f       Parse barcodes from Illumina headers.
chrom=f         For mapped sam files, make one file per chromosome (scaffold) using the rname.
header=f        Use the entire sequence header.
delimiter=      For prefix or suffix mode, specifying a delimiter will allow exact matches even if the length is variable.
                This allows demultiplexing based on names that are found without specifying a list of names.
                In suffix mode, for example, everything after the last delimiter will be used.
                Normally the delimiter will be used as a literal string (a Java regular expression); for example, ':' or 'HISEQ'.
                But there are some special delimiters which will be replaced by the symbol they name, 
                because they are reserved in some operating systems or cause other problems.
                These are provided for convenience due to possible OS conflicts:
                   space, tab, whitespace, pound, greaterthan, lessthan, equals,
                   colon, semicolon, bang, and, quote, singlequote
                These are provided because they interfere with Java regular expression syntax:
                   backslash, hat, dollar, dot, pipe, questionmark, star,
                   plus, openparen, closeparen, opensquare, opencurly
                In other words, to match '.', you should set 'delimiter=dot'.
substring=f     Names can be substrings of read headers.  Substring mode is
                slow if the list of names is large.  Requires a list of names.

Other Processing Parameters:
column=-1       If positive, split the header on a delimiter and match that column (1-based).
                For example, using this header:
                NB501886:61:HL3GMAFXX:1:11101:10717:1140 1:N:0:ACTGAGC+ATTAGAC
                You could demux by tile (11101) using 'delimiter=: column=5'
                Column is 1-based (first column is 1).
                If column is omitted when a delimiter is present, prefixmode
                will use the first substring, and suffixmode will use the last substring.
names=          List of strings (or files containing strings) to parse from read names.
                If the names are in text files, there should be one name per line.
                This is optional.  If a list of names is provided, files will only be created for those names.
                For example, 'prefixmode=t length=5' would create a file for every unique last 5 characters in read names,
                and every read would be written to one of those files.  But if there was addionally 'names=ABCDE,FGHIJ' 
                then at most 2 files would be created, and anything not matching those names would go to outu.
length=0        If positive, use a suffix or prefix of this length from read name instead of or in addition to the list of names.
                For example, you could create files based on the first 8 characters of read names.
hdist=0         Allow a hamming distance for demultiplexing barcodes.  This requires a list of names (barcodes).
replace=        Replace some characters in the output filenames.  For example, replace=+-
                would replace the + symbol in headers with the - symbol in filenames.  So you could
                match the name ACTGAGC+ATTAGAC in the header, but write to a file named ACTGAGC-ATTAGAC.

Buffering Parameters
streams=4       Allow at most this many active streams.  The actual number of open files
                will be 1 greater than this if outu is set, and doubled if output
                is paired and written in twin files instead of interleaved.
minreads=0      Don't create a file for fewer than this many reads; instead, send them to unknown.
                This option will incur additional memory usage.

Common parameters:
ow=t            (overwrite) Overwrites files that already exist.
zl=4            (ziplevel) Set compression level, 1 (low) to 9 (max).
int=auto        (interleaved) Determines whether INPUT file is considered interleaved.
qin=auto        ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto.
qout=auto       ASCII offset for output quality.  May be 33 (Sanger), 64 (Illumina), or auto (same as input).
                    

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                The max is typically 85% of physical memory.
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
	freeRam 3200m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

function demuxbyname() {
	local CMD="java $EA $EOOM $z $z2 -cp $CP jgi.DemuxByName2 $@"
	echo $CMD >&2
	eval $CMD
}

demuxbyname "$@"
