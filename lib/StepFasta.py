"""
Script running the fasta analysis step.
"""

import logging

import os
import FastaResFunc
import Init
from Logging import setup_logger

if __name__ == "__main__":
    # Init and run analysis steps
    snakemake = globals()["snakemake"]

    # Logging
    logger = logging.getLogger("main.fasta")
    setup_logger(logger, snakemake.log[0])

    config = snakemake.config
    config["queryName"] = str(snakemake.wildcards).split(":", 1)[0]
    config["output"] = str(snakemake.output)

    config["input"] = os.path.join(config["outdir"],config["queryName"]+"_accessions.txt")

    config["step"] = snakemake.rule

    parameters = Init.paramDef(config)

    fin = open(str(config["input"]), "r")
    lBlastres = list(map(str.strip, fin.readlines()))
    fin.close()

    # # Run step
    FastaResFunc.fastaCreation(parameters, lBlastres, parameters["output"])
