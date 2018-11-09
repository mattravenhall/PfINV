# PfINV
Detection of inversions in _P. falciparum_ from long read assemblies.

## Overview
This repo accompanies an upcoming publication describing the identification of inversions in long read _P. falciparum_ assemblies. It contains the core inversion detection script (src/), its support scripts (support/), and a number of scripts associated with the manuscript itself (manuscript/). The core discovery program is a command line Unix tool, built in bash and tested in both Linux and MacOS environments. It relies on `nucmer`, which can be automatically installed with `setup.sh`.

## Quick Start
```bash
# first time setup
bash setup.sh

# discovery with a list of sample assemblies (one path per line in samples.txt)
bash findINVs.sh -r=ref.fasta -s=samples.txt
```

## Method
Inversion detection relies upon `nucmer` alignment of a candidate fasta (the 'sample') against a reference fasta. Mis-alignments are then examined with `show-coords`, with those regions that align in reverse being identified as inversions. Aligned regions below a minimum size (default 300 bp) or whose relative bp coordinates (between the reference and sample assemblies) exceed a maximum distance (default 50,000 bp) will be discarded.

### Core Algorithm
1. `show-coords -HT ${SAMPLE}.delta`: Call `show-coords` on post-`nucmer` sample-on-reference alignment.
2. `awk -F'\t' '$3 > $4'`: Filter to aligned regions whose sample start location exceeds its sample end location (ie. it is in reverse alignment).
3. `sed 's,'"$REFname"'_,'"$SAMPLE"'_,g'`: Rename reference chromosome to sample's name for later matching.
4. `awk -F'\t' '$8 == $9'`: Exclude alignments in different chromosomes.
5. `sort -r -k5n`: Reverse sort alignments by the sample sequence length.  
6. `awk -v minSize="$MIN_SIZE" -F'\t' '$5 >= minSize'`: Remove alignments below or equal to a minimum size.
7. `awk -v maxDist="$MAX_DISTANCE" -F'\t' 'function abs(v) {return v < 0 ? -v : v} {if (abs($2 - $3) <= maxDist) {print $0}}'`: Remove alignments above a maximum relative distance (coords in sample vs reference).
8. `cut -f1-8`: Remove the final (now duplicate) column.
9. `> ${SAMPLE}.INV`: Write inversions to `sample`.INV.

## Output
Following a successful run `findINVs.sh` will output one .INV file per sample, in which each row represents an inversion. This file is formatted as:

StartRef | EndRef | StartSample | EndSample | LengthRef | LengthSample | Identity | SampleChrom
-------- | ------ | ----------- | --------- | --------- | ------------ | -------- | ----------
Reference start position (bp) | Reference end position (bp) | Sample start position (bp) | Sample end position (bp) | Length in reference | Length in sample | Percentage Identity | Chromosome |        

## Specific Scripts
### src/
- `setup.sh` Run to ensure that all required packages are installed. This will prompt to run as sudo.
- `findINVs.sh` Core detection pipeline, run as `findINVs.sh -h` to return the help page. Note that input fastas should be named &lt;sample>.fasta, and chromosome or contig names should be formatted as &lt;sample>_&lt;chromosome>. Input arguments include:
  - `-r --ref (reference file)` The path to your reference genome (.fasta).
  - `-s --samples (sample list)` The path to a textfile containing paths (one per line) to your sample fasta files.
  - `-m --minsize (min. size)` The minimum size for inversions, any below this cut-off will be discarded.
  - `-d --maxdist (max. distance)` The maximum relative distance in bp coordinates between your reference and candidate fastas.

### support/
- `grabRef.sh` Download the references used in the manuscript.
- `embl2fasta.py` Convert embl files to fasta format (for input to `findINVs.sh`).

### manuscript/
Scripts used to create the figures in the manuscript, free for use and adaptation by the wider community. Their filenames are indicative of their function.

## Raising Issues
If you find bugs, please raise them as an [issue on github](https://github.com/mattravenhall/PfINV/issues/new). For questions regarding the manuscript, email me at matt.ravenhall@lshtm.ac.uk.
