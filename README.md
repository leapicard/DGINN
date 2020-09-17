# DGINN: Detection of Genetic INNovations pipeline

DGINN is a pipeline dedicated to the detection of genetic innovations, starting from a nucleotidic sequence. 

It automatizes all the necessary preliminary steps for evolutionary analyses, including retrieval of homologs, 
assignment to orthology groups, codon alignment and reconstruction of gene phylogeny. 
Following the obtention of the alignements and corresponding phylogenies, three major genetic innovations are detected: 
duplication events, recombination events, and signatures of positive selection. 

DGINN was validated on nineteen primate genes with known evolutionary histories, and results can be consulted on BioRxiv 
(doi: https://doi.org/10.1101/2020.02.25.964155).
Results from the validation are available in the [corresponding repository](https://github.com/leapicard/DGINN_validation).
The version of DGINN used in the paper refers to [commit 5db0253](https://github.com/leapicard/DGINN/commit/5db02532408afcafad50a0b70dcf247ab4800492) 
and can be fetched through:
```{sh}
git init
git remote add origin https://github.com/leapicard/DGINN
git fetch --depth 1 origin 5db02532408afcafad50a0b70dcf247ab4800492
git checkout FETCH_HEAD
```
The docker is available for both the paper version and the current version of DGINN.

Any questions or suggestions about the program can be addressed to lea.picard [at] ens-lyon.fr, 
laurent.gueguen [at] univ-lyon1.fr or lucie.etienne [at] ens-lyon.fr.

# Overview

![alt text](https://github.com/leapicard/DGINN/blob/master/etc/pipeline_diagram.pdf)

# Installation

## 1/ Necessary dependencies and softwares

- Softwares and versions: [EMBOSS:6.6](http://en.bio-soft.net/format/emboss.html), [PhyML 3.0](https://github.com/stephaneguindon/phyml), [PRANK v.170427](http://wasabiapp.org/software/prank/prank_installation/), [Treerecs v1.0](https://gitlab.inria.fr/Phylophile/Treerecs), [HYPHY 2.3](http://www.hyphy.org/installation/), [Bio++ v.2](https://github.com/BioPP)
- Python (>3.5) and packages: Biopython, ete3, collections, logging, shlex, os, numpy, scipy, requests, pandas, statistics, time, re, argparse

## 2/ Docker

A [docker image](https://hub.docker.com/repository/docker/leapicard/dginn) is available to provide a way to use DGINN without requiring installation of any software except for [Docker](https://docs.docker.com/install/).

Get the docker:
```{sh}
docker pull leapicard/dginn
```

To download a specific version of the docker:
```{sh}
docker pull leapicard/dginn:paper
```

Use the docker:
```{sh}
docker run --rm -u $(id -u $USER):$(id -u $USER) -e HOME=. -v $PWD:$PWD -w $PWD leapicard/dginn -p <filename>
```
The command should be run as is, and should work on both Mac and Linux systems, provided the user belong to the 'docker' group (please refer to the [Docker Documentation](https://docs.docker.com/install/linux/linux-postinstall/) for help about setting the user as part of this group on Linux.)

All other arguments are passed exactly as if DGINN were run through the command line directly from the script (such as -p parameters.txt / see next section). However, one main difference is that all the files should be referred to by their name in the parameter file and be located within the working directory, while they can be referred by their path and be located in a different directory when running the script version.

# Usage

## 1/ Command line

```
DGINN, a pipeline for the Detection of Genetic Innovations.

optional arguments:
  -h, --help            show this help message and exit
  -dd, --debug          Enter verbose/debug mode

Mandatory parameters:
  -p <filename>, --params <filename>
                        Mandatory file with all the parameters necessary to
                        run the pipeline.

Optional parameters:
  -i <filename>, --infile <filename>
                        Path or list of paths (absolute or relative) to the
                        file(s) needed to start the pipeline (if indicated,
                        will take priority over the parameters file)
  -q <string>, --query <string>
                        Full identifier of the query in the format
                        SpeciesName_GeneName_GeneID (if indicated, will take
                        priority over the parameters file)
  -o <path>, --outdir <path>
                        Path to the output directory (if indicated, will take
                        priority over the parameters file)
  -host <filename>, --hostfile <filename>
                        Path to cluster hostfile if needed for mpi process
```

## 2/ Parameter file

DGINN uses a parameter file to pass all the necessary arguments for launching the pipeline.
Two example files are provided in the examples directory:
1. one performing steps 1-7 (see Overview) from the CDS of the gene of interest to the detection of recombination (parameters.txt)
2. one performing step 8 for the detection of positive selection (parameters_possel.txt)

This is the recommended usage for DGINN, so that analyses for positive selection can be parallelized over all alignments instead of doing them sequentially.

Please be aware that fasta sequence name **and** queryName must follow
the format speSpe_GENE_Id (ex: homSap_MX1_CCDS13673,
macMul_APOBEC3G_NM_001198693).

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

# Option for eliminating overly long sequences (default: cutoff(3))
# IQR or cutoff, factor can be put after in parenthesis
# cutoff will delete all sequences longer than (factor) times the median of the distribution
# IQR will delete all sequences longer than the third quartile plus (factor) times the InterQuartile Range
maxLen:

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

# Options for running PhyML
# Input the command in the same way you would to run PhyML yourself in the following manner phyml -i ALN [the rest of your options]
# For example, to run PhyML with a GTR model with 100 bootstraps, the option would be phymlOpt:phyml -i ALN -m GTR -b 100
# Please be aware that PhyML will run even if your options are wrong, but with its own default parameters
phymlOpt:

# Path to the species tree for the detection of duplication events and ortholog group assignment
# Species names must be formated as speSpe or speSpeSpe (ex: homSap, gorGorGor)
sptree:

# Option for the identification of duplication events (default: False)
duplication:

# Option for Long Branch separation (default: cutoff(50))
# IQR or cutoff, factor can be put after in parenthesis (ex: cutoff(50))
# EXPERIMENTAL
LBopt:

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

| Step              | Necessary file\(s\)                          | Format                     |
|-------------------|----------------------------------------------|----------------------------|
| blast             | CDS of the gene of interest                  | Fasta                      |
| accessions        | List of blast results                        | NCBI tabulated format (tsv)|
| fasta             | List of accession identifiers \(one/line\)   | Txt                        |
| orf               | mRNA sequences of orthologs                  | Fasta                      |
| alignment         | CDS sequences of orthologs                   | Fasta                      |
| tree              | \(codon\) alignment of orthologs             | Fasta                      |
| duplication       | \(codon\) alignment, gene tree               | Fasta, newick              |
| recombination     | \(codon\) alignment                          | Fasta                      |
| positiveSelection | codon alignment, gene tree                   | Fasta, gene tree           |

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

Users wishing to do the fastest check possible on their genes of interest are encouraged to run only BIO++ site models, 
as our validation results point to their providing the best compromise of solid results and shorter running times.

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

DGINN was validated on nineteen primate genes with known evolutionary histories, and results can be consulted on BioRxiv 
(doi: https://doi.org/10.1101/2020.02.25.964155).
Results from the validation are available in the [corresponding repository](https://github.com/leapicard/DGINN_validation).

# Utility scripts

## 1/ CCDSquery

In the etc folder, a script entitled CCDSquery.py is included. 
This script allows the user to download the CCDS sequences of human genes, by providing the properly formatted file obtained through HGNC.
This file should at least contain a column titled "Approved symbol" and another titled "CCDS accession".

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

Another script called parseResults.py can also be found in the etc folder. 

The input file is composed of two tab-separated columns: the first one indicates the full path to the directories containing the positive selection results (the directory containing the subdirectories busted, bpp_site, paml_site, etc.), the second one the full path to the alignments on which those analyses were performed.

Ex: /PATH/TO/GENENAME_sequences_filtered_longestORFs_mafft_mincov_prank_results_TIMESTAMP1/positive_selection_results_TIMESTAMP2 /PATH/TO/GENENAME_sequences_filtered_longestORFs_mafft_mincov_prank.best.fas

The script will output 3 different files: 
1. a summary of results (one gene per line)
2. the percentages of coverage at each position of each alignment (NB: it is advised not to modify this file in any capacity to ensure proper visualization of the results)
3. the likelihoods calculated by Bio++ (Bpp) and PAML codeml for each gene.

The different output files obtained with this script can be used to generate figures similar to those exposed in the DGINN paper through the [Shiny app](https://leapicard.shinyapps.io/DGINN-visualization/), which documentation can be found on the [corresponding repository](https://github.com/leapicard/DGINN-visualization).

```
python3 etc/parseResults.py -h
usage: etc/parseResults.py [-h] [-v] -in <filename>
                                      [-o <path/to/directory>]
                                      [-pr <value>]

This program outputs a summary of the results obtained through running DGINN
on a list of genes of interest.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         display /home/lea/Documents/DGINN/etc/parseResults.py
                        version number and exit

Mandatory input infos for running:
  -in <filename>, --inFile <filename>
                        List of all the directories containing the results
                        from DGINN analyses on different genes, and their
                        corresponding alignments.

Optional input infos (default values):
  -o <path/to/directory>, --outdir <path/to/directory>
                        folder for analysis results (path - by default output
                        file will be saved in the incoming directory)
  -pr <value>, --postrate <value>
                        folder for analysis results (path - by default output
                        file will be saved in the incoming directory)

```

