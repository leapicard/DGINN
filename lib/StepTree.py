"""
Script running the tree analysis step.
"""

import Init, os
import AnalysisFunc

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

    dAlTree = AnalysisFunc.runPhyML(parameters)

    os.rename(dAlTree, parameters["output"])
