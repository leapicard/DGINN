"""
Script running the ORF analysis step.
"""

import Init
import AnalysisFunc

if __name__ == "__main__":
    # Init and run analysis steps
    snakemake = globals()["snakemake"]

    config = snakemake.config
    config["queryName"] = str(snakemake.wildcards).split(":",1)[0]
    config["output"] = str(snakemake.output)
    cq = config["allquery"][config["queryName"]]
    if len(snakemake.input)>1 and cq!="void":
      config["input"] =  cq
    else:
      config["input"] = str(snakemake.input)
    config["step"] = snakemake.rule
    
    parameters = Init.paramDef(config)

    # Run step

    AnalysisFunc.getORFs(parameters)

