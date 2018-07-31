# PfINV
Detection of inversions in _P. falciparum_ from long read assemblies.

## Overview
This repo accompanies an [upcoming publication](MANUSCRIPT_LINK) describing the identification of inversions in long read _P. falciparum_ assemblies. It contains the core inversion detection script (src/), its support scripts (support/), and a number of scripts associated with the manuscript itself (manuscript/). The core discovery program is a command line Unix tool, built in bash and tested in both Linux and MacOS environments. It relies on `nucmer`, which can be automatically installed with `setup.sh`.

## Quick Run
```bash
# first time setup
bash setup.sh

# discovery with samples.txt
bash findINVs.sh -r=ref.fasta -s=samples.txt
```

## src/
- `setup.sh` Run to ensure that all required packages are installed. This will prompt to run as sudo.
- `findINVs.sh` Core detection pipeline, run as `findINVs.sh -h` to return the help page. Note that chromosome/contig names should be formatted as &lt;sample>_&lt;chromosome>, with fastas named &lt;sample>.fasta.

## support/
- `grabRef.sh` Download the references used for the [manuscript](MANUSCRIPT_LINK).
- `embl2fasta.py` Convert embl files to fasta format (for input to `findINVs.sh`).

## manuscript/
Scripts used to create the figures in the [manuscript](MANUSCRIPT_LINK), free for use and adaptation by the wider community. Their filenames are indicative of their function.

### Raising Issues
If you find bugs, please raise them as a [github issue](https://github.com/mattravenhall/PfINV/issues/new). For questions regarding the manuscript, email me at matt.ravenhall@lshtm.ac.uk.
