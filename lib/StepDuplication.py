"""
Script running the duplication analysis step.
"""

import logging
import os
import subprocess

import AnalysisFunc, TreeFunc
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

    align=config["input"].split()[0]
    
    config["step"] = snakemake.rule

    # else:
    parameters = Init.paramDef(config)

    # Run step

    parameters["sptree"] = TreeFunc.treeCheck(parameters.get("sptree",""), align, logger)
        
    dqueryaln = TreeFunc.splitTree(parameters)

    # output of the resulting sub-alignments querynames

    fout = open(config["output"], "w")

    if len(dqueryaln) >= 1:  # several sub alignments
        for query,faln in dqueryaln.items():
            fout.write(query + "\t"+ faln + "\n")

    else:
        fout.write(config["queryName"] + "\t"+ str(snakemake.input[0]) + "\n")
      
    fout.close()
