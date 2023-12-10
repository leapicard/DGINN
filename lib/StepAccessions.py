"""
Script running the accessions analysis step.
"""

import LoadFileFunc, Init
import ExtractFunc

if __name__ == "__main__":
    # Init and run analysis steps
    snakemake = globals()["snakemake"]

    config = snakemake.config
    config["output"] = str(snakemake.output)
    config["input"] = str(snakemake.input)
    config["queryName"] = str(snakemake.wildcards).split(":",1)[0]
    config["step"] = snakemake.rule
    
    parameters = Init.paramDef(config)

    # Run step
    query, lBlastRes = LoadFileFunc.accnEntry(parameters["input"])

    ExtractFunc.makeAccnsFile(lBlastRes, query, parameters["output"])

