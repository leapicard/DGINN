"""
Script running the positive selection analysis step.
"""

import Init
import PosSelFunc
import yaml, os, subprocess

if __name__ == "__main__":
    # Init and run analysis steps
    snakemake = globals()["snakemake"]

    config = snakemake.config
    config["queryName"] = str(snakemake.wildcards).split(":",1)[0]
    config["output"] = str(snakemake.output)
    config["step"] = snakemake.rule
    cq = config["allquery"][config["queryName"]]
    if len(snakemake.input)>1 and cq!="void" and config["step"] != "positive_selection":
      config["input"] =  cq
    else:
      config["input"] = str(snakemake.input)

    # Run step

    if config["input"].endswith("recombinations.txt") or config["input"].endswith("duplications.txt"):
        frec = open(config["input"],"r")
        lq=list(map(str.strip,frec.readlines()))
        frec.close()

        dpar={k:v for k,v in config.items() if k not in ["input","output","step"]}
        dpar["recombination"]=False
        dpar["duplication"]=False
        dpar["queryName"]=lq
        newconfig = "."+config["queryName"]+"config_sel.yaml"
        with open(newconfig,"w") as lout:
           yaml.dump(dpar,lout)

        subprocess.run(['snakemake',"--nolock","--cores=%d"%(max(1,len(lq))),"--configfile="+newconfig,"--until=positive_selection"])

        os.remove(newconfig)

        if len(lq)>1:  ## several subqueries, otherwise only queryName
          f=open(config["output"],"w")
          for quer in lq:
            fquer = open(config["outdir"] + "/" + quer + "_positive_selection.txt","r")
            for l in fquer:
              f.write(l)
            fquer.close()
            os.remove(config["outdir"] + "/" + quer + "_positive_selection.txt")
          f.close()

    else:
        parameters = Init.paramDef(config)
    
        outDir = PosSelFunc.pspAnalysis(parameters)
        aln = parameters["outdir"]+"/"+parameters["queryName"]+"_align.fasta"

        ## register resulting files in output
        f=open(config["output"],"w")
        f.write(outDir + "\t" + aln + "\n")
        f.close()


