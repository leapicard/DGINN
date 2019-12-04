# DGINN: Detection of Genetic INNovations pipeline

DGINN is a pipeline dedicated to the detection of genetic innovations, starting from a nucleotidic sequence. 

It automatizes all the necessary preliminary steps for evolutionary analyses, including retrieval of homologs, assignment to orthology groups, codon alignment and reconstruction of gene phylogeny.

Following the obtention of the alignements and corresponding phylogenies, three major genetic innovations: duplication events, recombination events, and signatures of positive selection.

Any questions or suggestions about the program can be addressed to lea.picard [at] ens-lyon.org

# Overview

![alt text](https://github.com/leapicard/DGINN/blob/master/etc/pipeline_diagram.png)

# Installation

## 1/ Necessary dependencies and softwares

- Softwares: [EMBOSS:6.6](http://en.bio-soft.net/format/emboss.html), [PhyML 3.0](https://github.com/stephaneguindon/phyml), [PRANK v.170427](http://wasabiapp.org/software/prank/prank_installation/), [Treerecs v1.0](https://gitlab.inria.fr/Phylophile/Treerecs), [HYPHY 2.3](http://www.hyphy.org/installation/), [Bio++ v.2](https://github.com/BioPP)
- Python (>3.5) and packages: Biopython, ete3, collections, logging, shlex, os, numpy, scipy, requests, pandas, statistics, time, re, argparse

## 2/ Docker

A docker image is available to provide a way to use DGINN without requiring installation of any software except for [Docker](https://docs.docker.com/install/).

Get the docker:
```{sh}
docker pull leapicard/dginn
```
Use the docker:
```{sh}
docker run -v /path/to/working/directory:/data leapicard/dginn -p parameters.txt
```
/path/to/working/directory should be the complete path to the directory where all files necessary to run the pipeline are located (parameter file, infile(s), species tree, etc.)

/data refers to the working directory within the docker and should not be changed.

All other arguments are passed exactly as if DGINN were run through the command line directly from the script (such as -p parameters.txt / see next section). However, one main difference is that all the files should be referred to by their name in the parameter file and be located within the working directory, while they can be referred by their path and be located in a different directory when running the script version.

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
                        Path or list of paths (absolute or relative) to the file(s) needed to start
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
Two example files are provided in the examples directory:
1. one performing steps 1-7 (cf Overview) from the CDS of the gene of interest to the detection of recombination (parameters.txt)
2. one performing step 8 for the detection of positive selection (parameters_possel.txt)

This is the recommended usage for DGINN, so that analyses for positive selection can be parallelized over all alignments instead of doing them sequentially.

Please be aware that fasta sequence names/queryName must follow the format speSpe_GENE_Id (ex: homSap_MX1_CCDS13673, macMul_APOBEC3G_NM_001198693).

```
# Path or list of paths (absolute or relative) to the files needed to start the pipeline
# Please refer to **3/ Entry steps** for necessary files
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

# Determines if Blast is performed against NCBI databases (default: True)
remote:

# NCBI API key to increase Blast speed, obtainable from the NCBI
APIKey:

# Path to the species tree for the detection of duplication events and ortholog group assignment
# Species names must be formated as speSpe or speSpeSpe (ex: homSap, gorGorGor)
sptree:

# Option for the identification of duplication events (default: False)
duplication:

# Minimum number of species for assignment to an ortholog group (default: 8)
nbspecies:

# Option for the detection of recombination events (default: False)
recombination:

# Option for the detection of positive selection (default: False)
positiveSelection:

# P-value for Hyphy methods (BUSTED/MEME) (Pond *et al.*, 2005) (default: 0.1)
hyphySeuil:

# Option for using the Hyphy method BUSTED (Murrel *et al.*, 2015) (default: False)
busted:

# Option for using the Hyphy method BUSTED (Murrel *et al.*, 2015) (default: False)
meme:

# Models to be computed by BIO++ (Gueguen *et al.*, 2013) and PAML (Yang, 2007)
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

| Step              | Necessary file\(s\)                          | Format                |
|-------------------|----------------------------------------------|-----------------------|
| blast             | CDS of the gene of interest                  | Fasta                 |
| accession         | List of blast results                        | NCBI tabulated format |
| fasta             | List of accession identifiers \(one/line\)   | Txt                   |
| orf               | mRNA sequences of orthologs                  | Fasta                 |
| alignment         | CDS sequences of orthologs                   | Fasta                 |
| tree              | \(codon\) alignment of orthologs             | Fasta                 |
| duplication       | \(codon\) alignment, gene tree, species tree | Fasta, newick, newick |
| recombination     | \(codon\) alignment                          | Fasta                 |
| positiveSelection | codon alignment, gene tree                   | Fasta, gene tree      |

File order must be respected and follow the one indicated in this table.


Though codon alignments are not technically necessary for the phyml, duplication and recombination steps, they are for positiveSelection.
Thus, starting at steps upstream of positiveSelection with non codon alignments will probably lead to failure at the positiveSelection step.

## 4/ Positive selection

DGINN includes different softwares to check for positive selection:
* BUSTED (Murrel et al., Molecular Biology and Evolution, 2015) from Hyphy
* MEME (Murrel et al., PLoS Genetics, 2012) from Hyphy
* PAML codeml (Yang, Molecular Biology and Evolution, 2007) for site models M0, M1, M2, M7 and M8
* BIO++ (Guéguen et al., Molecular Biology and Evolution, 2013) for site models M0, M1, M2, M7 and M8
* BIO++ for one-per-branch (OPB) model (similar to PAML codeml FreeRatio model) to test positive selection on branches

The first three methods are automatically parameterized in DGINN.

For BIO++, the parameter files can be automatically generated by DGINN, but the user can also provide their own parameter files if they wish to tweak the parameters further. The OPB option can also be used for different analyses using Bio++ as its results do not influence any subsequent step. Example parameter files for bppml and bppmixedlikelihoods (for site models) are provided in examples/, as well as a parameter file for running a one-per-branch model.

# Tutorial

## 1/ Example files

In the examples folder, two parameter files are provided.

NB: these files should be updated with the absolute paths to the files referred to instead of just their name when using DGINN through the command line and not through the docker.

```python3 DGINN.py -p parameters.txt```

Will launch DGINN steps 1-7 on ex_CCDS.fasta by :
* retrieving homologs of primate species in the NCBI *nr* database
* detecting duplications and assigning ortholog groups of at least 8 species based on ex_spTree.tree
* detecting recombination events

```python3 DGINN.py -p parameters_possel.txt```

Will launch DGINN step 8 on ex_aln.fasta and ex_genetree.tree by :
* looking for positive selection on the gene using BUSTED
* looking for sites under episodic positive selection using MEME
* looking for sites under positive selection using models M0-NS, M1-NS, M2-NS, M7-NS and M8-NS from BIO++
* looking for sites under positive selection using models M0, M1, M2, M7 and M8 from PAML codeml
* looking for branches under positive selection using BIO++

## 2/ Validation data

All results from our validation of DGINN will be made available shortly.

# Utility scripts

## 1/ CCDSquery

In the etc folder, a script entitled CCDSquery.py is included. 
This script allows the user to download the CCDS sequences of human genes, by providing the properly formatted file obtained through HGNC.
This file should at least contain a column titled "" and another titled "", as exemplified in examples/ex_CCDStable.txt.

```
python3 DGINN/etc/CCDSQuery.py -h
usage: DGINN/etc/CCDSQuery.py [-h] -in <filename>

This program get sequences' genes from HGNC Biomart.

optional arguments:
  -h, --help            show this help message and exit

Mandatory input infos for running:
  -in <filename>, --inFile <filename>
                        Table of HGNC approved symbols (one per line) and
                        corresponding CCDS accessions for the genes of
                        interest, obtained from HGNC Biomart
```

## 2/ Results extraction

Another script called parseResults.py can also be found in the etc folder. It compiles all the results from DGINN found in the directory passed as argument and outputs a summary of them.

```
python3 DGINN/etc/parseResults.py -h
usage: DGINN/etc/parseResults.py [-h] [-v] [-dd] -in <filename>
                                 [-o <path/to/directory>]

This program creates a summary of DGINN's results.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         display DGINN/etc/parseResults.py version number and
                        exit
  -dd, --debug          enter verbose/debug mode

Mandatory input infos for running:
  -in <filename>, --inDir <filename>
                        Path to directory of results

Optional input infos (default values):
  -o <path/to/directory>, --outdir <path/to/directory>
                        folder for analysis results (path - by default output
                        file will be saved in the incoming directory)

```

