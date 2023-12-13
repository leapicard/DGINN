import shutil, sys, os

# --- Path functions ---


def out_path(file, queryName = "{queryName}"):
    return expand("{outdir}/{queryName}_" + file, outdir=config["outdir"], queryName=queryName)


def log_path(file, queryName = "{queryName}"):
    return expand("{outdir}/logs/{queryName}_" + file, outdir=config["outdir"], queryName=queryName)


def data_path(file, queryName = "{queryName}"):
    return expand("{outdir}/data/{queryName}_" + file, outdir=config["outdir"], queryName=queryName)

          
####
# get queryName


queryName=config.get("queryName","")

if not queryName or len(queryName)==0:
   if "infile" in config and config["infile"]!="":
      if not "step" in config or config["step"]=="blast":
        name = config["infile"].split(",")[0].split("/")[-1]
        pn=name.rfind(".")
        if pn!=-1:
           queryName = name[:pn].strip()
        else:
           queryName = name.strip()
      else:
        try:
          queryName = os.path.split(config["infile"].split(",")[0])[-1].rsplit(".",1)[0].strip()
        except:
          print("Missing accession file " + config["infile"][0])
          sys.exit()        
   else:
     queryName="query"

          
# --- Default rule ---

up2ps = ("positiveSelection" in config and config["positiveSelection"])
with_dupli = ("duplication" in config and config["duplication"])
with_rec = ("recombination" in config and config["recombination"])

          
rule all:
    input:
        out_path("results.txt",queryName) if up2ps else
          (out_path("duplications.txt",queryName) if with_dupli else
            (out_path("recombinations.txt",queryName) if with_rec else
               (out_path("tree.dnd",queryName))))

               

# --- Step rules ---

## Needed starting rule for intermediate steps 
rule step:
    output:
        config["infile"],
          

rule blast:
    input:
        config["infile"],
    output:
        out_path("blastres.tsv"),
    log:
         log_path("01_blast.log"),
    script:
        "lib/StepBlast.py"
       
step = config.get("step","")

rule accessions:
    input:
        rules.blast.output if not step=="accessions" else config["infile"],
    output:
        out_path("accessions.txt"),
    log:
        log_path("02_accessions.log"),
    script:
        "lib/StepAccessions.py"


rule fasta:
    input:
       rules.accessions.output,
    output:
        out_path("sequences.fasta"),
    log:
        log_path("03_fasta.log"),
    script:
        "lib/StepFasta.py"


rule orf:
    input:
        rules.fasta.output if not step=="orf" else config["infile"],
    output:
        out_path("longestORFs.fasta"),
    log:
        log_path("04_orf.log"),
    script:
        "lib/StepORF.py"


          
rule alignment:
    input:
        rules.orf.output if not step=="alignment" else config["infile"],
    output:
        out_path("align.fasta"),
    log:
        log_path("05_alignment.log"),
    script:
        "lib/StepAlignment.py"


rule tree:
    input:
        out_path("align.fasta",queryName) if not step=="tree" else config["infile"],
    output:
        out_path("tree.dnd"),
    log:
        log_path("06_tree.log"),
    script:
        "lib/StepTree.py"


rule recombination:
    input:
        out_path("align.fasta") if not step=="recombation" else config["infile"],
    output:
        out_path("recombinations.txt"),
    log:
        log_path("08_recombination.log"),
    script:
        "lib/StepRecombination.py"


rule duplication:
    input:
        rules.recombination.output if with_rec else ([rules.alignment.output,rules.tree.output] if not step=="duplication" else config["infile"]),
    output:
        out_path("duplications.txt"),
    log:
        log_path("07_duplication.log"),
    script:
        "lib/StepDuplication.py"

rule positive_selection:
    input:
        rules.duplication.output if with_dupli else (rules.recombination.output if with_rec else ([rules.alignment.output,rules.tree.output] if not step=="positive_selection" else config["infile"])),
    output:
        out_path("positive_selection.txt"),
    log:
        log_path("08_positive_selection.log"),
    script:
        "lib/StepPositiveSelection.py"


rule analyse_ps:
    input:
        rules.positive_selection.output,
    output:
        out_path("results.txt"),
    log:
        log_path("09_results.log"),
    script:
        "lib/StepParseResults.py"
