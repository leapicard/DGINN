#!/bin/sh

mamba create -n dginn
mamba activate dginn
mamba install python=3.10
mamba install -c bioconda biopython mafft emboss ete3 prank phyml treerecs bpp-core bpp-seq bpp-phyl bpp-popgen hyphy paml snakemake
mamba install pandas