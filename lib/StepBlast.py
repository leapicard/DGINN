"""
Script running the blast analysis step.
"""
import logging
import os

import BlastFunc
import Init
from Logging import setup_logger

if __name__ == "__main__":
    # Init and run analysis steps
    snakemake = globals()["snakemake"]

    # Logging
    logger = logging.getLogger("main.blast")
    setup_logger(logger, snakemake.log[0])

    config = snakemake.config

    config["queryName"] = str(snakemake.wildcards).split(":", 1)[0]
    config["output"] = str(snakemake.output)

    config["input"] = os.path.join(config["outdir"],config["queryName"]+"_raw.fasta")

    config["step"] = snakemake.rule

    parameters = Init.paramDef(config)

    BlastFunc.blast(parameters, str(snakemake.output))
