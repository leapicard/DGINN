"""
Script running the ORF analysis step.
"""
import logging

import os
import AnalysisFunc
import Init
from Logging import setup_logger

if __name__ == "__main__":
    # Init and run analysis steps
    snakemake = globals()["snakemake"]

    # Logging
    logger = logging.getLogger("main.orf")
    setup_logger(logger, snakemake.log[0])

    config = snakemake.config
    config["queryName"] = str(snakemake.wildcards).split(":", 1)[0]
    config["output"] = str(snakemake.output)

    config["input"] = os.path.join(config["outdir"],config["queryName"]+"_sequences.fasta")

    config["step"] = snakemake.rule

    parameters = Init.paramDef(config)

    # Run step

    AnalysisFunc.getORFs(parameters)
