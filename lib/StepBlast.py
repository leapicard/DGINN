"""
Script running the blast analysis step.
"""

import Init
import BlastFunc

if __name__ == "__main__":
    # Init and run analysis steps
    snakemake = globals()["snakemake"]

    queryName = str(snakemake.wildcards).split(":",1)[0]
    
    config = snakemake.config
    config["output"] = str(snakemake.output)
    config["input"] = str(snakemake.input)
    config["queryName"] = queryName
    config["step"] = snakemake.rule
    
    parameters = Init.paramDef(config)

    BlastFunc.blast(parameters, str(snakemake.output))

