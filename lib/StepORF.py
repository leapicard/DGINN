"""
Script running the ORF analysis step.
"""
import logging

import AnalysisFunc
import Init
from Logging import setup_logger

if __name__ == "__main__":
    # Init and run analysis steps
    snakemake = globals()["snakemake"]

    # Logging
    logger = logging.getLogger("main.orf")
    setup_logger(logger, snakemake.log[0])

    config = snakemake.config
    config["queryName"] = str(snakemake.wildcards).split(":", 1)[0]
    config["output"] = str(snakemake.output)
    config["input"] = list(snakemake.input)

    cq = config["allquery"][config["queryName"]]
    if len(config["input"]) > 1 and cq != "void":
        config["input"] = cq
    else:
        config["input"] = str(snakemake.input)
    config["step"] = snakemake.rule

    parameters = Init.paramDef(config)

    # Run step

    AnalysisFunc.getORFs(parameters)
