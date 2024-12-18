# DGINN: Detection of Genetic INNovations pipeline

DGINN is a pipeline dedicated to the detection of genetic innovations, starting from a gene sequence.

It automatizes all the necessary preliminary steps for evolutionary
analyses, including retrieval of homologs, assignment to orthology
groups, codon alignment and reconstruction of gene phylogeny.

Once the alignements and corresponding phylogenies are obtained, three
major genetic innovations are detected: duplication events,
recombination events, and signatures of positive selection.

Eventually, the results of the whole process are summarized in a file `..._results.txt`.


DGINN was validated on nineteen primate genes with known evolutionary histories, and results can be consulted in the associated paper
(doi: https://doi.org/10.1093/nar/gkaa680).
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

![Diagram of DGINN main steps](https://github.com/leapicard/DGINN/blob/master/etc/pipeline_diagram.pdf)

# Installation

The pipeline is organized using snakemake.

## 1/ Necessary dependencies and softwares

-   Softwares and versions: [EMBOSS:6.6](http://en.bio-soft.net/format/emboss.html), [PRANK v.170427](http://wasabiapp.org/software/prank/prank_installation/), [MACSE V2.07](https://bio.tools/macse), [PhyML 3.0](https://github.com/stephaneguindon/phyml), [IQTREE v 2.6.6](http://www.iqtree.org), [Treerecs v1.0](https://gitlab.inria.fr/Phylophile/Treerecs), [HYPHY 2.3](http://www.hyphy.org/installation/), [Bio++ v.3](https://github.com/BioPP)
-   Python (>3.5) and packages: Biopython, ete3, collections, logging, shlex, os, numpy, scipy, requests, pandas, statistics, time, re, argparse
-   Snakemake

## 2/ Containers

The simplest way to use DGINN is through the use of a container, which
frees the user from the necessity of installing all of DGINN's
dependencies, and should make cross-platform usage possible (Linux/Mac
OS/Windows).

<!-- ### 2.1/ Conda -->

<!-- We provide a conda setting for DGINN, through file environment.yml. To -->
<!-- run locally, first create and activate a conda environment from -->
<!-- environment.yml: -->

<!-- ```shell -->
<!-- conda env create --file=environment.yml -->
<!-- conda activate dginn -->
<!-- ``` -->

<!-- then, in the working directory -->

<!-- ```shell -->
<!-- snakemake -s path_to_Snakefile --cores 1 --configfile=configuration_file -->
<!-- ``` -->

<!-- In case you want to allow more than one core for the analysis, set up -->
<!-- the "--cores" option accordingly. If the number is omitted (i.e., only -->
<!-- --cores is given), the number of used cores is determined as the -->
<!-- number of available CPU cores in the machine. -->


The user can use either of the images that we provide through
[Docker](https://docs.docker.com/install/) or
[Apptainer](https://apptainer.org/docs/user/latest/quick_start.html),
so the only software installation needed is the one for the chosen
container system.

Please be aware that, due to Docker compulsory root access, a Docker
container is not appropriate for usage in cluster environments, though
it is appropriate for cloud computing (tutorial to come) and local
usage. The Apptainer container should be usable in every environment.

#### a/ Docker

To use Docker, either you pull the image from github:

```{sh}
docker pull ghcr.io/lgueguen/dginn:master
```

or you to clone this repository first, and then build the Docker image
with:

```{sh}
docker build . -t dginn
```

After that, you will be able to run DGINN with:

```{sh}
docker run --rm -u $(id -u $USER) --network host --mount source=$(pwd),target=/opt/home,type=bind dginn --cores 1 --configfile config_example.yaml
```

The command should run on both Mac and Linux systems, provided the
user belongs to the 'docker' group (please refer to the [Docker
Documentation](https://docs.docker.com/install/linux/linux-postinstall/)
for help about setting the user as part of this group on Linux.)

We unfortunately cannot promise about the Docker container's usability
on Windows. In case the container doesn't work, we advise the user to
try the Apptainer container.

#### b/ Apptainer / Singularity

To use an Apptainer container you can pull it from github:

```{sh}
apptainer pull oras://ghcr.io/lgueguen/dginn:latest
```

To use the container, you can then run the following, for example with
"config_example.yaml" as your configuration file :

```{sh}
apptainer run --bind .snakemake:/opt/DGINN/.snakemake dginn_latest.sif --cores 1 --configfile config_example.yaml
```

<!-- If you want to run DGINN from another folder, you can specify the -->
<!-- path to the Snakefile file in the cloned repository: -->

<!-- ```{sh} -->
<!-- apptainer run --bind .snakemake:/home/mambauser/.snakemake /path/to/dginn.sif -s /path/to/Snakefile --cores 1 --configfile config_example.yaml -->

<!-- ``` -->

# Usage

## 1/ Configuration file

DGINN uses a parameter file in yaml syntax to pass all the necessary
arguments. Two example files are provided in the examples directory:

1. one performing steps 1-7 (see Overview) from the CDS of the gene of
   interest to the detection of recombination (parameters.yaml)

2. one performing step 8 for the detection of positive selection
   (parameters_possel.yaml) and summarizing the results in a file.

Please be aware that, when provided, fasta sequence name **and**
queryName must follow the format **speSpe_GENE_Id** for matching (ex:
homSap_MX1_CCDS13673, macMul_APOBEC3G_NM_001198693).

```
# Path or list of paths (absolute or relative) to the files needed to start the pipeline
infile:

# Output directory for all results
# Automatically created if not specified
outdir:

# Path to a file where progress of the pipeline will be logged
# Automatically created if not specified
logfile:

##################################
### STEP
##################################

# Step at which to enter the pipeline (default: blast)
# Please refer to 3/ Entry steps for names and necessary files
step:

##################################
### BLAST
##################################

# NCBI database on which the blast is to be performed (ex: nr)
# Future implementations will include the possibility to perform the search on local databases
blastdb:

# E-value for Blast (default: 10⁻⁴)
evalue:

# Coverage for Blast (default: 50)
mincov:

# Percentage of identity for Blast (default: 70)
percID:

#################################
### QUERY
#################################

# Option for eliminating overly long sequences (default: cutoff(3))
# IQR or cutoff, factor can be put after in parenthesis
# cutoff will delete all sequences longer than (factor) times the median of the distribution
# IQR will delete all sequences longer than the third quartile plus (factor) times the InterQuartile Range
maxLen:

# Can be used to limit the search on NCBI databases to certain set of species, to exclude others, etc.
# https://www.ncbi.nlm.nih.gov/books/NBK3837/#EntrezHelp.Entrez_Searching_Options
entryQuery:

# Identifier of the reference sequences (or list of identifiers of the
# same length as infile). Automatically filled is left empty
queryName:

# Determines if Blast is performed against NCBI databases (default: True)
remote:

# NCBI API key to increase Blast speed, obtainable from the NCBI
APIKey:

##################################################
###### ALIGNMENT
##################################################

# Choice of codon aligner: prank or macse (default)

aligner:

##################################################
###### TREE
##################################################

# Choice of tree builder: iqtree or phyml (default)

builder:


# Options for running PhyML
# Input the command in the same way you would to run PhyML yourself in the following manner phyml -i ALN [the rest of your options]
# For example, to run PhyML with a GTR model with 100 bootstraps, the option would be phymlOpt:phyml -i ALN -m GTR -b 100
# Please be aware that PhyML will run even if your options are wrong, but with its own default parameters
phymlOpt:


##################################################
###### ORTHOLOGS
##################################################

# Path to the species tree for the detection of duplication events and ortholog group assignment
# Species names must be formated as speSpe or speSpeSpe (ex: homSap, gorGorGor)
sptree:

# Option for the identification of duplication events (default: False)
duplication:

###############################################
##### CLEANING
###############################################

# Option for Long Branch separation (default: cutoff(50))
# IQR or cutoff, factor can be put after in parenthesis (ex: cutoff(50))
# EXPERIMENTAL
LBopt:

# Minimum number of species for assignment to an ortholog group (default: 8)
nbspecies:

##################################################
###### RECOMBINATION
##################################################

# Option for the detection of recombination events (default: False)
recombination:

##################################################
###### POSITIVE SELECTION
##################################################

# Option for the detection of positive selection (default: False)
positiveSelection:

# P-value for Hyphy methods (BUSTED/MEME) (Pond *et al.*, 2005) (default: 0.1)
hyphySeuil:

# Option for using the Hyphy method BUSTED (Murrel *et al.*, 2015) (default: False)
busted:

# Option for using the Hyphy method BUSTED (Murrel *et al.*, 2015) (default: False)
meme:

# Models to be computed by BIO++ (Gueguen *et al.*, 2013) and/or PAML (Yang, 2007)
# Implemented models: M0, M1, M2, M7, M8, M8a, DFP07, DFP07_0
# Must be comma separated (ex: M0,M1,M2)
#
# Rate distribution are either Constant ou Gamma(n=4)
# Default is Gamma with Bio++ and Constant with PAML, and explicit
# rate distribution are available through "_C" or "_G" suffixes to model
# names (ex: M0_C, M0_G)

models: M0, M1, M2

# Option for using paml for the detection of sites under positive selection (default: False)
paml:

# Option for using BIO++ for the detection of sites under positive selection
# If True, parameter file will be automatically generated
# Can be used to indicate the path to a BIO++ parameter file
bppml:

# Same as previously, but for extracting results from the results computed from bppml
mixedlikelihood:

# Option for using BIO++ for the detection of branches under positive selection
# If True, parameter file will be automatically generated
# Positive selection on each is assessed through LRT M2 vs M1 model in bio++.
# Parameters different from omega are shared between all branches.
opb:
```

## 2/ Entry steps

| Step              | Necessary file\(s\)                        | Format                      |
| ----------------- | ------------------------------------------ | --------------------------- |
| blast             | CDS of the gene of interest                | Fasta                       |
| accessions        | List of blast results                      | NCBI tabulated format (tsv) |
| fasta             | List of accession identifiers \(one/line\) | Txt                         |
| orf               | mRNA sequences of orthologs                | Fasta                       |
| alignment         | CDS sequences of orthologs                 | Fasta                       |
| tree              | \(codon\) alignment of orthologs           | Fasta                       |
| duplication       | \(codon\) alignment, gene tree             | Fasta, newick               |
| recombination     | \(codon\) alignment                        | Fasta                       |
| positiveSelection | codon alignment, gene tree                 | Fasta, gene tree            |

File order must be respected and follow the one indicated in this table.

Though alignments are not technically necessary in codons for the
tree, duplication and recombination steps, they are for
positiveSelection.

## 3/ Positive selection

DGINN includes different softwares to check for positive selection:

-   BUSTED (Murrel et al., Molecular Biology and Evolution, 2015) from Hyphy
-   MEME (Murrel et al., PLoS Genetics, 2012) from Hyphy
-   PAML codeml (Yang, Molecular Biology and Evolution, 2007) for site models M0, M1, M2, M7, M8a and M8
-   BIO++ (Guéguen et al., Molecular Biology and Evolution, 2013) for site models M0, M1, M2, M7 M8a, M8, M10 and DFP_07
-   BIO++ for one-per-branch (OPB) model (similar to PAML codeml FreeRatio model) to test positive selection on branches

The two first methods are automatically parameterized in DGINN.

For BIO++, the parameter files can be automatically generated by
DGINN, but the user can also provide their own parameter files if they
wish to tweak the parameters further. 
<!-- The OPB option can also be used -->
<!-- for different analyses using Bio++ as its results do not influence any -->
<!-- subsequent step.  -->
Example parameter files for bppml and
bppmixedlikelihoods (for site models) are provided in examples, as
well as a parameter file for running a one-per-branch model.

Users wishing to do the fastest check possible on their genes of
interest are encouraged to run only BIO++ site models (with Constant
rate distribution), as our validation results point to their providing
the best compromise of solid results and shorter running times.

## 4/ Output

# Tutorial

## 1/ Example files

In the examples folder, two parameter files are provided.

NB: these files should be updated with the paths to the files referred
to instead of just their name when using DGINN through the command
line and not through the docker.

In the next command, it is provided that the user is in `examples` directory.

```
snakemake -s ../Snakefile --cores 1 --configfile=parameters.yaml
```

Will launch DGINN steps 1-7 on ex_CCDS.fasta by :

-   retrieving homologs of primate species in the NCBI _nr_ database
-   detecting duplications and assigning ortholog groups of at least 8 species based on ex_spTree.tree
-   detecting recombination events

```
snakemake -s ../Snakefile --cores 1 --configfile=parameters_poosel.yaml
```

Will launch DGINN step 8 on ex_aln.fasta and ex_genetree.tree by :

-   looking for positive selection on the gene using BUSTED
-   looking for sites under episodic positive selection using MEME
-   looking for sites under positive selection using models M0-NS, M1-NS, M2-NS, M7-NS, M8a-NS and M8-NS from BIO++
-   looking for sites under positive selection using models M0, M1, M2, M7, M8a and M8 from PAML codeml
-   looking for branches under positive selection using BIO++

## 2/ Validation data

DGINN was validated on nineteen primate genes with known evolutionary histories, and results can be consulted on BioRxiv
(doi: https://doi.org/10.1101/2020.02.25.964155).
Results from the validation are available in the [corresponding repository](https://github.com/leapicard/DGINN_validation).

# Utility scripts

Several utility scripts upstream and downstream of DGINN in the etc folder,

## 1/ Upstream

CCDSquery.py allows the user to download the CCDS sequences of human
genes, by providing the properly formatted file obtained through HGNC.
This file should at least contain a column titled "Approved symbol"
and another titled "CCDS accession".

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

## 2/ Downstream

These programs are useful if the last rule `analyse_ps` has not run
(which should not occur).

### 1/ recup_to_parse

recup_to_parse produces a file that will be used by parseResult (see
below) to parse the output of DGINN analyses.

recup_to_parse is to be run in the directory from which all DGINN
analyses have been run, and where the output files
"tag_DGINN_date.log" are written. The script scans all those file,
keeping the most recent successful one for each tag, in case there
have been several analyses.

Message outputs whether analyses seem to have worked well (+tag) or
positive selection but no clean exit (~tag), or not at all (-tag).
Any gene with + or ~ sign is written in the parsing file.

```
python3 recup_to_parse.py [-o outfile]

where:

-o outfile: output file that will parsed by parseResults.py

```

### 2/ parseResult

parseResult parses a file describing where results of DGINN are to be
found for several genes, and output a summary of the analyses in a
single file.

The input file is composed of two tab-separated columns: the first one
indicates the full path to the directories containing the positive
selection results (the directory containing the subdirectories busted,
bpp_site, paml_site, etc.), the second one the full path to the
alignments on which those analyses were performed.

Ex: /PATH/TO/GENENAME_sequences_filtered_longestORFs_mafft_mincov_prank_results_TIMESTAMP1/positive_selection_results_TIMESTAMP2 /PATH/TO/GENENAME_sequences_filtered_longestORFs_mafft_mincov_prank.best.fas

The script will output 2 different files:

1. a summary of results (one gene per line)
2. the percentages of coverage at each position of each alignment (NB: it is advised not to modify this file in any capacity to ensure proper visualization of the results)
 <!--- 3. the likelihoods calculated by Bio++ (Bpp) and PAML codeml for each gene. --->

The different output files obtained with this script can be used to generate figures similar to those exposed in the DGINN paper through the [Shiny app](https://leapicard.shinyapps.io/DGINN-visualization/), which documentation can be found on the [corresponding repository](https://github.com/leapicard/DGINN-visualization).

```
python3 DGINN/etc/parseResults.py -h
usage: DGINN/etc/parseResults.py [-h] [-v] -in <filename> [-o <path/to/directory>] [-pr <value>] [-pm <value>]

This program outputs a summary of the results obtained through running DGINN on a list of genes of interest.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         display DGINN/etc/parseResults.py version number and exit

Mandatory input infos for running:
  -in <filename>, --inFile <filename>
                        List of all the directories containing the results from DGINN analyses on different genes, and their corresponding alignments.

Optional input infos (default values):
  -o <path/to/directory>, --outdir <path/to/directory>
                        folder for analysis results (path - by default output file will be saved in the incoming directory)
  -pr <value>, --postrate <value>
                        Threshold posterior probability of omega>1 to admit positive selected sites.
  -pm <value>, --pvmeme <value>
                        Maximum p-value of PS site significance for MEME method.cd
```

# Citation

In case of usage of DGINN, please cite:
Lea Picard, Quentin Ganivet, Omran Allatif, Andrea Cimarelli, Laurent Guéguen, Lucie Etienne, DGINN, an automated and highly-flexible pipeline for the detection of genetic innovations on protein-coding genes, Nucleic Acids Research, Volume 48, Issue 18, 09 October 2020, Page e103, https://doi.org/10.1093/nar/gkaa680
