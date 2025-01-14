"""
Script running the recombination analysis step.
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
    logger = logging.getLogger("main.recombination")
    setup_logger(logger, snakemake.log[0])

    config = snakemake.config
    config["queryName"] = str(snakemake.wildcards).split(":", 1)[0]
    config["output"] = str(snakemake.output)
 
    config["input"] = os.path.join(config["outdir"],config["queryName"]+"_align.fasta")

    parameters = Init.paramDef(config)

    # Run step (inactivated for now)

    lbp = AnalysisFunc.runPhymlMulti(parameters)
    
    gardRes = os.path.join(config["outdir"],config["queryName"] + "_bp.txt")

    frec = open(gardRes, "w")
    frec.write("breakpoints\n")
    frec.write(str(lbp)+"\n")
    frec.close()

    ## lQuer is the list of new queryNames
    ## lAln is the list of new alignments
    [lQuer, lAln] = AnalysisFunc.parseGard(gardRes, parameters)

    ## register resulting [queryName, alignment]s in output
    f = open(config["output"], "w")
    for i in range(len(lQuer)):
        f.write(lQuer[i] + "\t" + lAln[i] + "\n")
        
    f.close()
