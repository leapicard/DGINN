import shutil, sys, os
from pathlib import Path
          
# --- Path functions ---
def out_path(file, queryName = "{queryName}"):
    return expand("{outdir}/{queryName}" + file, outdir=config["outdir"], queryName=queryName)


def log_path(file, queryName = "{queryName}"):
    return expand("{outdir}/logs/{queryName}" + file, outdir=config["outdir"], queryName=queryName)


def data_path(file, queryName = "{queryName}"):
    return expand("{outdir}/data/"+os.path.basename("{}".format(queryName)) + file, outdir=config["outdir"], queryName=queryName)


#### get step
          
step = config.get("step","blast")
if not step:
          step="blast"
          
## get Snakefile path

snakefilepath=workflow.source_path("Snakefile")
pfile=snakefilepath.rfind(os.sep+"file"+os.sep)
snakefilepath=snakefilepath[pfile+4+len(os.sep):]
          
####
# match infile & queryName

### for duplication, either list of queryNames, or couple of infile [align, tree]                    
if not "infile" in config or not config["infile"]:
    config["infile"]=["void"]

if type(config["infile"]) == str:
    config["infile"]=[config["infile"]]

## use absolute paths
config["outdir"]=os.path.abspath(config["outdir"])
if not os.path.exists(config["outdir"]):
          os.mkdir(config["outdir"])

if "sptree" in config and config["sptree"]:
   config["sptree"]=os.path.abspath(config["sptree"])

config["infile"]=[",".join([os.path.abspath(f2.strip()) if f2!="void" else "void" for f2 in f.split(",")]) for f in config["infile"]]
          
### in case of no queryNames, build them from infiles
if not "queryName" in config or not config["queryName"] or len(config["queryName"])==0:
  if config["infile"]!=["void"]:
    config["queryName"] = list(map(lambda x:os.path.split(x.split(",")[0])[-1].rsplit(".",1)[0].strip(), config["infile"]))
  else:
    print("Missing infile and queryName for duplication step")
    sys.exit(0)
elif type(config["queryName"]) == str:
    config["queryName"]=[config["queryName"]]


# Check everything is fine and match queryName:infile
if config["infile"]!=["void"]:
    if len(config["queryName"])!=len(config["infile"]):
       print("lengths of queryName & infile do not match.")
       sys.exit(0)


## adjust infiles to correct filenames and set up step
dstep={
          "blast":["_raw.fasta","_blastres.tsv"],
          "accessions":["_blastres.tsv","_accessions.txt"],
          "fasta":["_accessions.txt","_sequences.fasta"],
          "orf":["_sequences.fasta","_orf.fasta"],
          "alignment":["_orf.fasta","_align.fasta"],
          "tree":["_align.fasta","_tree.dnd"],
          "recombination":["_align.fasta","_recombination.txt"]
       }          

if step in dstep:
  for i in range(len(config["queryName"])):
    instep = dstep[step][0]
    shutil.copyfile(config["infile"][i],out_path(instep,queryName=config["queryName"][i])[0])
    config["infile"][i]=out_path(instep,queryName=config["queryName"][i])[0]


## specifically for positive selection step, infiles should be couples (alignment, tree)

if step == "positiveSelection" or step == "duplication":
  for i in range(len(config["queryName"])):
    infi = config["infile"][i].split(",")
    if len(infi)!=2:
           print(step + " step needs couples of infiles \"alignment, tree\", not " + config["infile"][i])
           sys.exit(0)
    outali = out_path(dstep["alignment"][1],queryName=config["queryName"][i])[0]
    shutil.copyfile(infi[0],outali)
    outtree = out_path(dstep["tree"][1],queryName=config["queryName"][i])[0]
    shutil.copyfile(infi[1],outtree)
    config["infile"][i]=",".join([outali,outtree])
          
### Go to snakefile directory
os.chdir(snakefilepath[:snakefilepath.rfind(os.sep)])

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


######################################################
#### Request
######################################################

## Function to check if a rule should be run: if file with outsuffix
## does not exists, ask for insuffix file
          
def check_exists(wildcards, outsuffix, insuffix):
    outfile = out_path(outsuffix,queryName=wildcards)[0]
    infile = out_path(insuffix,queryName=wildcards)[0]
    if not os.path.exists(outfile):
          return infile
    if os.path.getsize(outfile)==0:
          Path(infile).touch()
          return infile
    return ""
          
rule blast:
    input:
        out_path(dstep["blast"][0])
    output:
        out_path(dstep["blast"][1]),
    log:
        log_path("_01_blast.log"),
    script:
        "lib/StepBlast.py"
       
checkpoint accessions:
    input:
        lambda wc : check_exists(wc, "_accessions.txt", "_blastres.tsv"),
    output:
        out_path("_accessions.txt"),
    log:
        log_path("_02_accessions.log"),
    script:
        "lib/StepAccessions.py"


rule fasta:
    input:
        rules.accessions.output
    output:
        out_path("_sequences.fasta"),
    log:
        log_path("_03_fasta.log"),
    script:
        "lib/StepFasta.py"


checkpoint orf:
    input:
        lambda wc : check_exists(wc, "_orf.fasta", "_sequences.fasta"),
    output:
        out_path("_orf.fasta"),
    log:
        log_path("_04_orf.log"),
    script:
        "lib/StepORF.py"


######################################################
#### Alignment
######################################################

checkpoint alignment:
    input:
        lambda wc : check_exists(wc, "_align.fasta", "_orf.fasta"),
    output:
        out_path("_align.fasta"),
    log:
        log_path("_05_alignment.log"),
    script:
        "lib/StepAlignment.py"


######################################################
#### Tree
######################################################
          
checkpoint tree:
    input:
        lambda wc : check_exists(wc, "_tree.dnd", "_align.fasta"),
    output:
        out_path("_tree.dnd"),
    log:
        log_path("_06_tree.log"),
    script:
        "lib/StepTree.py"

          

######################################################
#### Recombination
######################################################

def get_list_from_ck(ck_out):
     faln = open(ck_out,"r")
     lqueryName = []
     for l in faln:
         lqueryName.append(l.split()[0].strip())
     faln.close()
     return lqueryName
          
def make_recombination(wildcards, suffix):
     ck_out = checkpoints.recombination_ck.get(**wildcards).output[0]
     return out_path(suffix,queryName=get_list_from_ck(ck_out))

          
checkpoint recombination_ck:
    input:
        out_path("_align.fasta"),
    output:
        out_path("_recombinations_ck.txt"),
    log:
        log_path("_08_recombination_ck.log"),
    script:
        "lib/StepRecombination.py"

rule recombination:
    input:     
       lambda wc: make_recombination(wc, "_align.fasta.txt")
    output:
       out_path("_recombinations.txt"),
    log:
       log_path("_08_recombination.log"),
    run:
       for nf in input:
          os.system("cat " + nf + " >> " + str(output))

                 
######################################################
#### Duplication
######################################################
                              
def make_duplication(wildcards, suffix = "_tree.dnd"):
     ck_out = checkpoints.duplication_ck.get(**wildcards).output[0]
     return out_path(suffix,queryName=get_list_from_ck(ck_out))

def make_rec_dup(wildcards, suffix):
     ck_out = checkpoints.recombination_ck.get(**wildcards).output[0]
     ls = get_list_from_ck(ck_out)
     ls2 = []
     for qN in ls:
        wildcards.queryName = qN
        ck2 = checkpoints.duplication_ck.get(**wildcards).output[0]
        ls2 += get_list_from_ck(ck2)
     return out_path(suffix,queryName=ls2)
          
checkpoint duplication_ck:
    input:
        [rules.alignment.output,rules.tree.output]
    output:
        out_path("_duplications_ck.txt"),
    log:
        log_path("_07_duplication_ck.log"),
    script:
        "lib/StepDuplication.py"


rule duplication_rec:
    input:
      make_duplication
    output:
      out_path("_duplications_rec.txt"),
    log:
      log_path("_07_duplication_rec.log"),
    run:
      for nf in input:
          os.system("cat " + nf + " >> " + str(output))


rule duplication:
    input:
       lambda wc: make_recombination(wc, "_duplications_rec.txt") if with_rec else make_duplication(wc)
    output:
       out_path("_duplications.txt"),
    log:
      log_path("_07_duplication.log"),
    run:
      for nf in input:
             os.system("cat " + nf + " >> " + str(output))


######################################################
#### Positive selection
######################################################

rule positive_selection_dup:
    input:
        [rules.alignment.output, rules.tree.output],
    output:
        out_path("_positive_selection_dup.txt"),
    log:
        log_path("_08_positive_selection_dup.log"),
    script:
        "lib/StepPositiveSelection.py"

rule positive_selection:
    input:
       ((lambda wc: make_rec_dup(wc, "_positive_selection_dup.txt")) if with_rec else (lambda wc: make_duplication(wc, "_positive_selection_dup.txt"))) if with_dupli else ((lambda wc: make_recombination(wc, "_positive_selection_dup.txt")) if with_rec else (lambda wc: out_path("_positive_selection_dup.txt"))),
    output:
        out_path("_positive_selection.txt"),
    log:
        log_path("_08_positive_selection.log"),
    run:
      for nf in input:
          os.system("cat " + nf + " >> " + str(output))



##################################################
### Analyse
##################################################
          
rule analyse_ps:
    input:
        rules.positive_selection.output,
    output:
        out_path("_results.txt"),
    log:
        log_path("_09_results.log"),
    script:
        "lib/StepParseResults.py"
