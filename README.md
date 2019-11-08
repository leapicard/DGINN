# DGINN: Detection of Genetic INNovations pipeline

IN CONSTRUCTION

DGINN is a pipeline dedicated to the detection of genetic innovations, starting from a nucleotidic sequence. 
It automatizes all the necessary preliminary steps for evolutionary analyses, including retrieval of homologs, assignment to orthology groups, codon alignment and reconstruction of gene phylogeny.
Following the obtention of the alignements and corresponding phylogenies, three major genetic innovations: duplication events, recombination events, and signatures of positive selection.
Any questions or suggestions about the program can be addressed to lea.picard [at] ens-lyon.org

# Installation

## 1/ Necessary dependencies and softwares

- Softwares: [Blast+ 2.6](https://www.ncbi.nlm.nih.gov/books/NBK279671/), [EMBOSS:6.6](http://en.bio-soft.net/format/emboss.html), [PhyML 3:3.1](https://github.com/stephaneguindon/phyml), [PRANK v.150803](http://wasabiapp.org/software/prank/prank_installation/), [Treerecs v1.0](https://gitlab.inria.fr/Phylophile/Treerecs), [HYPHY 2.220161014beta](http://www.hyphy.org/installation/), [PAML 4.9](http://abacus.gene.ucl.ac.uk/software/paml.html), [Bio++ v.2](https://github.com/BioPP)
- Python (>3.5) and packages: Biopython, ete3, collections, logging, shlex, os, numpy, scipy, requests, pandas, statistics, time, re, argparse

## 2/ Docker

A docker image will soon be available (link) to provide a way to use DGINN without requiring installation of any software except for [Docker](https://docs.docker.com/install/).

```{sh}
docker pull link/tocome
```

# Usage

## 1/ Command line

```
optional arguments:
  -h, --help            show this help message and exit
  -dd, --debug          Enter verbose/debug mode

Mandatory parameters:
  -p <filename>, --params <filename>
                        Mandatory file with all the parameters necessary to
                        run the pipeline.

Optional parameters:
  -i <filename>, --infile <filename>
                        Path or list of paths to the file(s) needed to start
                        the pipeline (if indicated, will take priority over
                        the parameters file)
  -q <string>, --query <string>
                        Full identifier of the query in the format
                        SpeciesName_GeneName_GeneID (if indicated, will take
                        priority over the parameters file)
  -host <filename>, --hostfile <filename>
                        Path to cluster hostfile if needed for mpi process
```

## 2/ Parameter file

DGINN uses a parameter file to pass all the necessary arguments for lauching the pipeline.
Two example files are provided in the main directory:
one performing steps 1-7 (cf pipeline workflow) from the CDS of the gene of interest to the detection of recombination,
and one performing step 8 for the detection of positive selection. 
This is the recommended usage for DGINN, so that analyses for positive selection can be parallelized over all alignments instead of doing them sequentially.
```
# Path or list of paths to the files needed to start the pipeline
# Please refer to 3/ Entry steps to know which files are necessary depending on the step at which DGINN will start
infile:

# NCBI database on which the blast is to be performed (ex: nr)
# Future implementations will include the possibility to perform the search on local databases
blastdb:

# Output directory for all results
# Automatically created if not specified
outdir:

# Path to a file where progress of the pipeline will be logged
# Automatically created if not specified
logfile:

# E-value for Blast (default: 10⁻⁴)
evalue:

# Coverage for Blast (default: 50)
mincov:

# Percentage of identity for Blast (default: 70)
percID:

# Can be used to limit the search on NCBI databases to certain set of species, to exclude others, etc.
# https://www.ncbi.nlm.nih.gov/books/NBK3837/#EntrezHelp.Entrez_Searching_Options
entryQuery:

# Step at which to enter the pipeline (default: blast)
# Please refer to 3/ Entry steps for names and necessary files
step:

# Identifier of the reference sequence for steps outside of blast and positiveSelection
queryName:

# Determines if Blast is performed against local or NCBI databases (default: True)
remote:

# NCBI API key to increase Blast speed, obtainable from the NCBI
APIKey:

# Path to the species tree for the detection of duplication events and ortholog group assignment
# Species names must be formated as speSpe or speSpeSpe (ex: homSap, gorGorGor)
sptree:

# Option for the identification of duplication events (default: False)
treerecs:

# Minimum number of species for a duplication to be considered as delimiting a group of orthologs (default: 8)
nbspecies:

# Option for the detection of recombination events (default: False)
gard:

# Option for the detection of positive selection (default: False)
positiveSelection:

# P-value for Hyphy methods (BUSTED/MEME) (Pond et al., 2005) (default: 0.1)
hyphySeuil:

# Option for using the Hyphy method BUSTED (Murrel et al., 2015) (default: False)
busted:

# Option for using the Hyphy method BUSTED (Murrel et al., 2015) (default: False)
meme:

# Models to be computed by BIO++ (Gueguen et al., 2013) and PAML (Yang, 2007)
# Implemented models: M0, M1, M2, M7, M8
# Must be comma separated (ex: M0,M1,M2)
models:

# Option for using BIO++ for the detection of sites under positive selection
# If True, parameter file will be automatically generated
# Can be used to indicate the path to a BIO++ parameter file
bppml:

# Same as previously, but for extracting results from the results computed from bppml
mixedlikelihood:

# Option for using BIO++ for the detection of branches under positive selection
# If True, parameter file will be automatically generated
# Can be used to indicate the path to a BIO++ parameter file
opb:
```

## 3/ Entry steps

## 4/ Positive selection

# Tutorial

## 1/ Example files

example sh script
example parameter file

## 2/ Validation data

link to validation results repo

# Related content

## 1/ CCDSquery

link to CCDSquery repo

## 2/ Results extraction

link to Results extraction repo


