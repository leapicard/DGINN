from snakemake.utils import min_version, validate
import logging

##### set minimum snakemake version #####
#min_version("6.15.1")

container: "docker://laugueguen/dginn"

##### setup report #####

configfile:"config/configfile.json"
config_file = "config/configfile.json"

validate(config, "config/config.schema.yaml")

#definition du dossier de stockage des résultats

parameters = config["parameters"]
data = config["data"]
resDir = data["o"]

filename = parameters["infile"].split("data/")[1].split(".")[0]
treename = parameters["treefile"].split("data/")[1].split(".")[0] if (parameters["step"]=="option") else None

logging.basicConfig(filename=f"results/{filename}.log",
                            level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s')


if parameters["step"]=="blast" and not parameters["positiveSelection"] and not parameters["duplication"] and not parameters["recombination"] :
    
    rule all:
      input: resDir+"tree_"+filename+"_filtered2.phylip_phyml_tree.txt"

    include: "rules/blast.smk"
    include: "rules/alignment.smk"
    include: "rules/phylogneticTree.smk"
  
  
elif parameters["step"]=="blast" and not parameters["positiveSelection"] and parameters["duplication"] and not parameters["recombination"] :
    
    rule all:
      input: resDir+"dup"+filename+"_recs2.nhx"
    
    include: "rules/blast.smk"
    include: "rules/alignment.smk"
    include: "rules/phylogneticTree.smk"
    include: "rules/duplication.smk"


elif parameters["step"]=="blast" and parameters["positiveSelection"] and not parameters["duplication"] and not parameters["recombination"] :
    rule all:
        input: resDir+"possel_"+filename+"_files_list.txt"

    include: "rules/blast.smk"
    include: "rules/alignment.smk"
    include: "rules/phylogneticTree.smk"
    include: "rules/positiveSelection.smk"


elif parameters["step"]=="blast" and parameters["positiveSelection"] and parameters["duplication"] and not parameters["recombination"] :
    rule all:
        input: resDir+"possel_"+filename+"_files_list.txt"
    
    include: "rules/blast.smk"
    include: "rules/alignment.smk"
    include: "rules/phylogneticTree.smk"
    include: "rules/duplication.smk"
    include: "rules/positiveSelection.smk"
  
elif parameters["step"]=="option" and parameters["duplication"] and parameters["positiveSelection"] and not parameters["recombination"]:
    rule all:
      input: resDir+"possel_"+filename+"_files_list.txt"
    
    include: "rules/duplication.smk"
    include: "rules/positiveSelection.smk"

elif parameters["step"]=="option" and parameters["duplication"] and not parameters["positiveSelection"] and not parameters["recombination"]:
    rule all:
      input: resDir+"dup_"+filename+"_recs2.nhx"
    
    include: "rules/duplication.smk"

elif parameters["step"]=="option" and parameters["positiveSelection"] and not parameters["duplication"] and not parameters["recombination"]:
    rule all:
      input: resDir+"possel_"+filename+"_files_list.txt"

    include: "rules/positiveSelection.smk"
  
elif parameters["step"]=="option" and not parameters["positiveSelection"] and not parameters["duplication"] and parameters["recombination"]:
    rule all:
      input: resDir+"recomb_"+filename+"_files_list.txt"

    include: "rules/recombination.smk"

elif parameters["step"]=="blast" and parameters["duplication"] and parameters["positiveSelection"] and parameters["recombination"]:

    rule all:
      input: resDir+"possel_"+filename+"_files_list.txt"
    
    include: "rules/blast.smk"
    include: "rules/mapping2.smk"
    include: "rules/phylogneticTree.smk"
    include: "rules/recombination.smk"
    include: "rules/duplication.smk"
    include: "rules/positiveSelection.smk"



    



   