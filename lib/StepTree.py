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

    builder = config.get("builder","phyml")

    config["input"] = os.path.join(config["outdir"],config["queryName"]+"_align.fasta")
    parameters = Init.paramDef(config)

    # Run step
    lbuilder=["phyml","iqtree"]

    treeOk=False
    while not treeOk:
      if builder == "phyml":
        logger.info("Running PhyML")
        dAltree = AnalysisFunc.runPhyML(parameters)
      elif builder == "iqtree":
        logger.info("Running IQtree")
        dAltree = AnalysisFunc.runIqTree(parameters)
      else:
        logger.info("Unknown tree builder: " + builder)
        break
   
      if not os.path.exists(dAltree) or os.path.getsize(dAltree)==0:
        logger.info(builder + " failed to build tree.")
        lbuilder = [b for b in lbuilder if b!=builder]
        if lbuilder==[]:
          break
        builder = lbuilder[0]
      else:
        treeOk=True
        break
   
    if not treeOk:
      raise Exception("Failed tree construction.")
   
    os.rename(dAltree, config["output"])
