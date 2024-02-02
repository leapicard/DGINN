"""
Script running the duplication analysis step.
"""

import logging
import os
import subprocess

import AnalysisFunc
import Init
import yaml
from Logging import setup_logger

if __name__ == "__main__":
    # Init and run analysis steps
    snakemake = globals()["snakemake"]

    # Logging
    logger = logging.getLogger("main.duplication")
    setup_logger(logger, snakemake.log[0])

    config = snakemake.config
    config["queryName"] = str(snakemake.wildcards).split(":", 1)[0]
    config["output"] = str(snakemake.output)
    config["input"] = str(snakemake.input)

    config["step"] = snakemake.rule

    # else:
    parameters = Init.paramDef(config)

    # Run step

    check_params = {p: parameters[p] for p in ["nbspecies", "LBopt"]}

    lquery = AnalysisFunc.checkTree(parameters)

    ## rerun snakemake if needed

    fout = open(config["output"], "w")

    if len(lquery) >= 1:  # several sub alignments
        for query in lquery:
            fout.write(query + "\t"+ config["outdir"] + "/" + query + "_longest.fasta" + "\n")

    else:
        fout.write(config["queryName"] + "\t"+ config["outdir"] + "/" + config["queryName"] + "_longest.fasta" + "\n")
      
    fout.close()
