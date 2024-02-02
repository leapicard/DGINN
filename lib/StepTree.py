"""
Script running the tree analysis step.
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
    logger = logging.getLogger("main.tree")
    setup_logger(logger, snakemake.log[0])

    config = snakemake.config

    config = snakemake.config
    config["queryName"] = str(snakemake.wildcards).split(":", 1)[0]
    config["output"] = str(snakemake.output)

    config["input"] = os.path.join(config["outdir"],config["queryName"]+"_align.fasta")
    parameters = Init.paramDef(config)

    # Run step

    dAlTree = AnalysisFunc.runPhyML(parameters)

    os.rename(dAlTree, config["output"])
