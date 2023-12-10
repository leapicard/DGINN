"""
Script running the recombination analysis step.
"""

import Init, subprocess
import AnalysisFunc
import yaml, os

if __name__ == "__main__":
    # Init and run analysis steps
    snakemake = globals()["snakemake"]

    config = snakemake.config
    config["output"] = str(snakemake.output)
    config["input"] = str(snakemake.input)
    config["queryName"] = str(snakemake.wildcards).split(":",1)[0]
    config["step"] = snakemake.rule
    
    parameters = Init.paramDef(config)

    # Run step (inactivated for now)

    gardRes = config["outdir"]+"/"+config["queryName"]+"_bp.txt"
    frec = open(gardRes,"w")
    frec.write("breakpoints\n")
    frec.write("[]\n")
    frec.close()
    
    lAln = AnalysisFunc.parseGard(gardRes, parameters)

    ## rerun snakemake if needed
    if len(lAln)>1: # several sub 
      dpar={k:v for k,v in config.items() if k not in ["infile","output","step"]}
      dpar["infile"]=""
      dpar["queryName"]=lAln
      dpar["recombination"]=False
      newconfig = "."+config["queryName"]+"config_dupl.yaml"
      with open(newconfig,"w") as lout:
        yaml.dump(dpar,lout)

      subprocess.run(['snakemake',"--cores=%d"%(len(lAln)),"--nolock","--configfile=" + newconfig,"--until=alignment"])
      
      os.remove(newconfig)

    ## register resulting files in output
    f=open(str(snakemake.output),"w")
    for aln in lAln:
      f.write(aln + "\n")
    f.close()
