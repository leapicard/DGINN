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
All
```
infile:
blastdb:
outdir:
logfile:
evalue:
mincov:
percID:
entryQuery:
step:
queryName:
remote:
APIKey:
sptree:
treerecs:
nbspecies:
gard:
positiveSelection:
hyphySeuil:
busted:
meme:
models:
bppml:
mixedlikelihood:
opb:
gnh:
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


