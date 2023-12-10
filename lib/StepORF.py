"""
Script running the ORF analysis step.
"""

import Init
import AnalysisFunc

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

    # Run step

    AnalysisFunc.getORFs(parameters)

