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


    config["queryName"] = str(snakemake.wildcards).split(":", 1)[0]
    config["output"] = str(snakemake.output)
    config["step"] = snakemake.rule

    builder = config.get("builder","phyml")

    config["input"] = os.path.join(config["outdir"],config["queryName"]+"_align.fasta")
    parameters = Init.paramDef(config)

    # Run step

    if builder == "phyml":
      dAltree = AnalysisFunc.runPhyML(parameters)
    elif builder == "iqtree":
      dAltree = AnalysisFunc.runIqTree(parameters)
    else:
      print("Unknown tree builder: " + builder)
      
    os.rename(dAltree, config["output"])
