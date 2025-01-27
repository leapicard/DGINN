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
    lbuilder=["phyml","iqtree"]

    flag=0

    while flag<=1:
      if builder == "phyml":
        dAltree = AnalysisFunc.runPhyML(parameters)
      elif builder == "iqtree":
        dAltree = ""#AnalysisFunc.runIqTree(parameters)
      else:
        logger.info("Unknown tree builder: " + builder)
        break
        
      if not os.path.exists(dAltree) or os.path.getsize(dAltree)==0:
        logger.info(builder + " failed to build tree.")
        lbuilder = [b for b in lbuilder if b!=builder]
        if lbuilder==[]:
          flag = 0
          break
        builder = lbuilder[0]
        flag+=1
      else:
        break

    if flag==0:
      raise Exception("Failed tree construction.")
    
    os.rename(dAltree, config["output"])
