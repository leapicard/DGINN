"""
Script running the fasta analysis step.
"""

import Init
import FastaResFunc

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
    
    fin=open(str(snakemake.input),"r")
    lBlastres = list(map(str.strip,fin.readlines()))
    fin.close()
    
    # # Run step
    FastaResFunc.fastaCreation(parameters, lBlastres, parameters["output"])

