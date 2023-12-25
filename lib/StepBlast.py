"""
Script running the blast analysis step.
"""

import Init
import BlastFunc

if __name__ == "__main__":
    # Init and run analysis steps
    snakemake = globals()["snakemake"]

    config = snakemake.config

    config["queryName"] = str(snakemake.wildcards).split(":",1)[0]
    config["output"] = str(snakemake.output)
    config["input"] = list(snakemake.input)

    cq = config["allquery"][config["queryName"]]
    if len(config["input"])>1 and cq!="void":  # if input is a list, take the correct element in it
      config["input"] =  cq
    else:
      config["input"] = str(snakemake.input)
      
    config["step"] = snakemake.rule

    parameters = Init.paramDef(config)

    BlastFunc.blast(parameters, str(snakemake.output))

