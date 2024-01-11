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
    config["queryName"] = str(snakemake.wildcards).split(":",1)[0]
    config["output"] = str(snakemake.output)
    config["input"] = str(snakemake.input)

    config["step"] = snakemake.rule

    if config["input"].endswith("recombinations.txt"):
        frec = open(config["input"],"r")
        lq=[]
        laln=[]
        for l in frec.readlines():
          ll=l.split()
          lq.append(ll[0].strip())
          laln.append(ll[1].strip())

        frec.close()

        dpar={k:v for k,v in config.items() if k not in ["input","output","step","infile"]}
        dpar["recombination"]=False
        dpar["queryName"]=lq
        dpar["infile"]=laln
        dpar["step"]="alignment"

        newconfig = "."+config["queryName"]+"_config_dupl.yaml"
        with open(newconfig,"w") as lout:
           yaml.dump(dpar,lout)

        subprocess.run(['snakemake',"--nolock","--cores=%d"%(max(len(lq),1)),"--configfile="+newconfig,"--until=duplication","--reason"])
        os.remove(newconfig)

        if len(lq)>1:  # several files, otherwise only first file
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

        fout=open(config["output"],"w")
    
        if len(lquery)>=1: # several sub alignments
          dpar={k:v for k,v in config.items() if k not in ["input","output","step"]}
          dpar["queryName"]=list(lquery.keys())
          dpar["infile"]=list(lquery.values())
          dpar["recombination"]=False
          dpar["step"]="alignment"
          
          newconfig = "."+config["queryName"]+"config_dupl.yaml"
          with open(newconfig,"w") as lout:
            yaml.dump(dpar,lout)
            
          subprocess.run(['snakemake',"--nolock","--cores=%d"%(max(1,len(lquery))),"--configfile="+ newconfig ,"--until=tree"])

          os.remove(newconfig)

          for query in lquery:
            fout.write(query + "\t"+ config["outdir"] + "/" + query + "_align.fasta" + "\n")
    
        else:
          fout.write(config["queryName"] + "\t"+ config["outdir"] + "/" + config["queryName"] + "_align.fasta" + "\n")
          
        fout.close()
    
