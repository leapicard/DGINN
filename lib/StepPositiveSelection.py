"""
Script running the positive selection analysis step.
"""
import logging
import os, subprocess

import Init
import PosSelFunc
import yaml
from Logging import setup_logger

if __name__ == "__main__":
    # Init and run analysis steps
    snakemake = globals()["snakemake"]

    # Logging
    logger = logging.getLogger("main.positiveSelection")
    setup_logger(logger, snakemake.log[0])

    config = snakemake.config
    config["queryName"] = str(snakemake.wildcards).split(":", 1)[0]
    config["output"] = str(snakemake.output)

    config["step"] = snakemake.rule
    aln = os.path.join(config["outdir"], config["queryName"] + "_align.fasta")

    # Run step
    parameters = Init.paramDef(config)

    outDir = PosSelFunc.pspAnalysis(parameters)

    ## register resulting files in output
    f = open(config["output"], "w")
    f.write(outDir + "\t" + aln + "\n")
    f.close()
