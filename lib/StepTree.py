"""
Script running the tree analysis step.
"""
import logging
import os, shutil

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

    config["queryName"] = str(snakemake.wildcards).split(":", 1)[0]
    config["output"] = str(snakemake.output)
    config["step"] = snakemake.rule

    if "builder" in config and config["builder"]:
        builder = [config["builder"]]
    else:
        builder = ["phyml","iqtree"]

    config["input"] = os.path.join(config["outdir"],config["queryName"]+"_align.fasta")
    parameters = Init.paramDef(config)

    # Run step

    lAltree = [] # List of (llist, file name)
    if "phyml" in builder:
        lAltree.append(AnalysisFunc.runPhyML(parameters))
        logger.info("Log-lik: " + str(lAltree[-1][0]))
    if "iqtree" in builder:
        lAltree.append(AnalysisFunc.runIqTree(parameters))
        logger.info("Log-lik: " + str(lAltree[-1][0]))

    ### get maximum
    if len(lAltree)==0:
      raise Exception("Failed tree construction.")

    imax = max(range(len(lAltree)), key = lambda x:lAltree[x][0])
    logger.info("Get tree from best builder: " + builder[imax])
    dAltree = lAltree[imax][1]
    
    os.rename(dAltree, config["output"])
