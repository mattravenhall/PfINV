#!/bin/bash

# Inversion Detection w/ nucmer

# Requires:
# - NUCMER (v3.1+)
# - show-coords
###

# Inputs:
# Reference fasta: Reference contig names must be named 'REF_<chrnum>', sample fasta contigs must be named '<SAMPLE>_<chrnum>'
# Samples: Formatted as <SAMPLE>_<CHR_NUM>.fa
# MIN_SIZE: Removes tiny inversions (default 300 bp)
# MAX_DISTANCE: Resolves for relative positions

# Returns
# <SAMPLE>.INV file, as 

# The Gist
# Perform whole chromosome alignment with nucmer, then isolate those regions within a certain buffer zone that are reversed

# Default values
REF='./reference.fa'
SAMPLES='./samples.txt'
MIN_SIZE=300		# Remove small inversions
MAX_DISTANCE=50000	# Min distance for reference and sample locations

# Help Info
help='Usage: findINVs.sh <--argument=value>\n
	-r --ref\tPath to reference fasta.\t(default: '${REF}')\n
	-s --samples\tPath to samples list.\t\t(default: '${SAMPLES}')\n
	-m --minsize\tMinimum size of inversion.\t(default: '${MIN_SIZE}')\n
	-d --maxdist\tMaximum relative distance.\t(default: '${MAX_DISTANCE}')'

# Parse arguments
if [ $# -eq 0 ]; then
	echo -e ${help}
	exit
fi 
for i in "$@"; do
	case $i in
		-h=*--help=*)
		echo -e ${help}
		exit
		;;
		-r=*|--ref=*)
		REF="${i#*=}"
		OLDIFS=$IFS
		IFS='.' read -ra ARRY <<< "$REF"
		REFname="${ARRY[0]}"
		IFS=$OLDIFS
		shift
		;;
		-s=*|--samples=*)
		readarray SAMPLES < ${i#*=}
		shift
		;;
		-m=*|--minsize=*)
		MIN_SIZE="${i#*=}"
		shift
		;;
		-d=*--maxdist=*)
		shift
		;;
		*)
		# Read as help
		echo -e ${help}
		exit
		;;
	esac
done

echo '=findINVs v1.0='
echo 'Reference: '${REF}
echo 'Samples: '${SAMPLES[@]}
echo 'Minimum Size: '${MIN_SIZE}
echo 'Maximum Distance: '${MAX_DISTANCE}

for SAMPLE in ${SAMPLES[@]}; do
	echo -e '\n=Aligning '${SAMPLE}'...'
	nucmer -p ${SAMPLE} ${REF} ./${SAMPLE}.fa*

	echo -e '\n=Isolating top inversions...'
	show-coords -HT ${SAMPLE}.delta | awk -F'\t' '$3 > $4' | sed 's,'"$REFname"'_,'"$SAMPLE"'_,g' | awk -F'\t' '$8 == $9' | sort -r -k5n | awk -v minSize="$MIN_SIZE" -F'\t' '$5 >= minSize' | awk -v maxDist="$MAX_DISTANCE" -F'\t' 'function abs(v) {return v < 0 ? -v : v} {if (abs($2 - $3) <= maxDist) {print $0}}' | cut -f1-8 > ${SAMPLE}.INV
	cat <(echo -e 'StartRef\tEndRef\tStartSample\tEndSample\tLengthRef\tLengthSample\tIdentity\tSampleChrom') ${SAMPLE}.INV | tee ${SAMPLE}.INV
done
