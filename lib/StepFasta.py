"""
Script running the fasta analysis step.
"""

import Init
import FastaResFunc

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
    
    fin=open(str(config["input"]),"r")
    lBlastres = list(map(str.strip,fin.readlines()))
    fin.close()
    
    # # Run step
    FastaResFunc.fastaCreation(parameters, lBlastres, parameters["output"])
