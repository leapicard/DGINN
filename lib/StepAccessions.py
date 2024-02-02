"""
Script running the accessions analysis step.
"""

import logging

import ExtractFunc
import Init
import LoadFileFunc
from Logging import setup_logger

if __name__ == "__main__":
    # Init and run analysis steps
    snakemake = globals()["snakemake"]

    # Logging
    logger = logging.getLogger("main.accessions")
    setup_logger(logger, snakemake.log[0])

    config = snakemake.config
    config["queryName"] = str(snakemake.wildcards).split(":", 1)[0]
    config["output"] = str(snakemake.output)

    config["input"] = str(snakemake.input)
    
    config["step"] = snakemake.rule

    parameters = Init.paramDef(config)

    # Run step
    query, lBlastRes = LoadFileFunc.accnEntry(parameters["input"])

    ExtractFunc.makeAccnsFile(lBlastRes, query, parameters["output"])
