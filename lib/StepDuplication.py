"""
Script running the duplication analysis step.
"""

import Init, os
import AnalysisFunc
import yaml, subprocess

if __name__ == "__main__":
    # Init and run analysis steps
    snakemake = globals()["snakemake"]

    config = snakemake.config
    config["output"] = str(snakemake.output)
    config["input"] = str(snakemake.input)
    config["queryName"] = str(snakemake.wildcards).split(":",1)[0]
    config["step"] = snakemake.rule

    if config["input"].endswith("recombinations.txt"):
        frec = open(config["input"],"r")
        lq=list(map(str.strip,frec.readlines()))
        frec.close()

        dpar={k:v for k,v in config.items() if k not in ["input","output","step"]}
        dpar["recombination"]=False
        dpar["queryName"]=lq

        newconfig = "."+config["queryName"]+"config_dupl.yaml"
        with open(newconfig,"w") as lout:
           yaml.dump(dpar,lout)

        subprocess.run(['snakemake',"--nolock","--cores=%d"%(len(lq)),"--configfile="+newconfig,"--until=duplication"])

        os.remove(newconfig)

        if len(lq)>1:  # several queries, otherwise only queryName
          f=open(str(snakemake.output),"w")
          for quer in lq:
            fquer = open(config["outdir"] + "/" + quer + "_duplications.txt","r")
            for l in fquer:
              f.write(l)
            fquer.close()
          os.remove(config["outdir"] + "/" + quer + "_duplications.txt")
          f.close()
        
    else:
        parameters = Init.paramDef(config)
    
        # Run step
    
        check_params = {p: parameters[p] for p in ["nbspecies", "LBopt"]}
    
        lquery = AnalysisFunc.checkTree(parameters)
        
        ## rerun snakemake if needed

        fout=open(str(snakemake.output),"w")
    
        if len(lquery)>1: # several sub alignments
          dpar={k:v for k,v in config.items() if k not in ["input","output","step"]}
          dpar["queryName"]=lquery
          dpar["recombination"]=False
          newconfig = "."+config["queryName"]+"config_dupl.yaml"
          with open(newconfig,"w") as lout:
            yaml.dump(dpar,lout)
            
          subprocess.run(['snakemake',"--nolock","--cores=%d"%(len(lquery)),"--configfile="+ newconfig ,"--until=tree"])

          os.remove(newconfig)

          for query in lquery:
            fout.write(query + "\n")
    
        else:
          fout.write(config["queryName"] + "\n")
          
        fout.close()
    
