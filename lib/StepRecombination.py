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
    config["queryName"] = str(snakemake.wildcards).split(":",1)[0]
    config["output"] = str(snakemake.output)
    cq = config["allquery"][config["queryName"]]
    if len(snakemake.input)>1 and cq!="void":
      config["input"] =  cq
    else:
      config["input"] = str(snakemake.input)
    config["step"] = snakemake.rule
    
    parameters = Init.paramDef(config)

    # Run step (inactivated for now)

    gardRes = config["outdir"]+"/"+config["queryName"]+"_bp.txt"
    frec = open(gardRes,"w")
    frec.write("breakpoints\n")
    frec.write("[500]\n")
    frec.close()

    ## lQuer is the list of new queryNames 
    lQuer = AnalysisFunc.parseGard(gardRes, parameters)
    
    ## rerun snakemake if needed
    if len(lQuer)>1: # several sub 
      dpar={k:v for k,v in config.items() if k not in ["output","step","queryName","infile"]}
      dpar["queryName"]=lQuer
      dpar["recombination"]=False
      newconfig = "."+config["queryName"]+"_config_rec.yaml"
      with open(newconfig,"w") as lout:
        yaml.dump(dpar,lout)

      subprocess.run(['snakemake',"--cores=%d"%(max(1,len(lQuer))),"--nolock","--configfile=" + newconfig,"--until=alignment"])
      
      os.remove(newconfig)

    ## register resulting queryNames in output
    f=open(config["output"],"w")
    for i in range(len(lQuer)):
      f.write(lQuer[i] + "\n")
      
    f.close()
