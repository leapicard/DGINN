import shutil, sys, os

# --- Path functions ---


def out_path(file, queryName = "{queryName}"):
    return expand("{outdir}/{queryName}" + file, outdir=config["outdir"], queryName=queryName)


def log_path(file, queryName = "{queryName}"):
    return expand("{outdir}/logs/{queryName}" + file, outdir=config["outdir"], queryName=queryName)


def data_path(file, queryName = "{queryName}"):
    return expand("{outdir}/data/"+os.path.basename("{}".format(queryName)) + file, outdir=config["outdir"], queryName=queryName)

          
####
# get queryName

if not "infile" in config or not config["infile"]:
    config["infile"]="void"

if type(config["infile"]) == str:
    config["infile"]=[config["infile"]]

if not "queryName" in config or not config["queryName"] or len(config["queryName"])==0:
    config["queryName"] = list(map(lambda x:os.path.split(x)[-1].rsplit(".",1)[0].strip(), config.get("infile","void")))
elif type(config["queryName"]) == str:
    config["queryName"]=[config["queryName"]]

# Check everything is fine (with the exception of positive_selection
# step, where two files (alignment & tree) are needed

if config["infile"]!=["void"] and config.get("step","")!="positive_selection":
  if len(config["queryName"])!=len(config["infile"]):
     print("lengths of queryName & infile do not match")
     sys.exit(0)
  else:
     config["allquery"]={config["queryName"][i]:config["infile"][i] for i in range(len(config["queryName"]))}
elif config["infile"]==["void"]:
     config["allquery"]={config["queryName"][i]:"void" for i in range(len(config["queryName"]))}
else: # In case of positive_selection step
     config["allquery"]={queryName:config["infile"]}
                    
# --- Default rule ---

up2ps = ("positiveSelection" in config and config["positiveSelection"])
with_dupli = ("duplication" in config and config["duplication"])
with_rec = ("recombination" in config and config["recombination"])

          
rule all:
    input:
        out_path("_results.txt",config["queryName"]) if up2ps else
          (out_path("_duplications.txt",config["queryName"]) if with_dupli else
            (out_path("_recombinations.txt",config["queryName"]) if with_rec else
               (out_path("_tree.dnd",config["queryName"]))))

               
# --- Step rules ---

## Needed starting rule for intermediate steps 
rule step:
    output:
     touch(expand("{filename}", filename=config["infile"]))

rule blast:
    input:
        rules.step.output,
    output:
        out_path("_blastres.tsv"),
    log:
        log_path("_01_blast.log"),
    script:
        "lib/StepBlast.py"
       
step = config.get("step","")

rule accessions:
    input:
        rules.blast.output if not step=="accessions" else rules.step.output,
    output:
        out_path("_accessions.txt"),
    log:
        log_path("_02_accessions.log"),
    script:
        "lib/StepAccessions.py"


rule fasta:
    input:
       rules.accessions.output,
    output:
        out_path("_sequences.fasta"),
    log:
        log_path("_03_fasta.log"),
    script:
        "lib/StepFasta.py"


rule orf:
    input:
        rules.fasta.output if not step=="orf" else rules.step.output,
    output:
        out_path("_longestORFs.fasta"),
    log:
        log_path("_04_orf.log"),
    script:
        "lib/StepORF.py"


          
rule alignment:
    input:
        rules.orf.output if not step=="alignment" else rules.step.output,
    output:
        out_path("_align.fasta"),
    log:
        log_path("_05_alignment.log"),
    script:
        "lib/StepAlignment.py"


rule tree:
    input:
        out_path("_align.fasta") if not step=="tree" else rules.step.output,
    output:
        out_path("_tree.dnd"),
    log:
        log_path("_06_tree.log"),
    script:
        "lib/StepTree.py"


rule recombination:
    input:
        out_path("_align.fasta") if not step=="recombination" else rules.step.output,
    output:
        out_path("_recombinations.txt"),
    log:
        log_path("_08_recombination.log"),
    script:
        "lib/StepRecombination.py"

rule duplication:
    input:
        rules.recombination.output if with_rec else ([rules.alignment.output,rules.tree.output] if not step=="duplication" else rules.step.output),
    output:
        out_path("_duplications.txt"),
    log:
        log_path("_07_duplication.log"),
    script:
        "lib/StepDuplication.py"

rule positive_selection:
    input:
        rules.duplication.output if with_dupli else (rules.recombination.output if with_rec else ([rules.alignment.output,rules.tree.output] if not step=="positive_selection" else rules.step.output)),
    output:
        out_path("_positive_selection.txt"),
    log:
        log_path("_08_positive_selection.log"),
    script:
        "lib/StepPositiveSelection.py"


rule analyse_ps:
    input:
        rules.positive_selection.output,
    output:
        out_path("_results.txt"),
    log:
        log_path("_09_results.log"),
    script:
        "lib/StepParseResults.py"
