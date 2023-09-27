# --- Main config file ---


# --- Path functions ---


def out_path(file):
    return expand("{outdir}/" + file, outdir=config["outdir"])


def log_path(file):
    return expand("{outdir}/logs/" + file, outdir=config["outdir"])


def data_path(file):
    return expand("{outdir}/data/" + file, outdir=config["outdir"])


# --- Default rule ---


rule all:
    input:
        data_path("data_blast.pkl"),
        data_path("data_accessions.pkl"),
        data_path("data_fasta.pkl"),
        data_path("data_orf.pkl"),
        data_path("data_alignment.pkl"),
        data_path("data_tree.pkl"),
        data_path("data_duplication.pkl"),
        data_path("dAlTree_duplication.pkl"),


# --- Step rules ---


rule blast:
    input:
        config["infile"],
    output:
        data_path("data_blast.pkl"),
        out_path("blastres.tsv"),
    log:
        log_path("01_blast.log"),
    script:
        "lib/StepBlast.py"


rule accessions:
    input:
        data_path("data_blast.pkl"),
        out_path("blastres.tsv"),
    output:
        data_path("data_accessions.pkl"),
        out_path("accessions.txt"),
    log:
        log_path("02_accessions.log"),
    script:
        "lib/StepAccessions.py"


rule fasta:
    input:
        data_path("data_accessions.pkl"),
        out_path("accessions.txt"),
    output:
        data_path("data_fasta.pkl"),
        out_path("sequences.fasta"),
    log:
        log_path("03_fasta.log"),
    script:
        "lib/StepFasta.py"


rule orf:
    input:
        data_path("data_fasta.pkl"),
        out_path("sequences.fasta"),
    output:
        data_path("data_orf.pkl"),
        out_path("allORFs.fasta"),
        out_path("longestORFs.fasta"),
    log:
        log_path("04_orf.log"),
    script:
        "lib/StepORF.py"


rule alignment:
    input:
        data_path("data_orf.pkl"),
        out_path("allORFs.fasta"),
        out_path("longestORFs.fasta"),
    output:
        data_path("data_alignment.pkl"),
        out_path("mafft.fasta"),
        out_path("clustiso.fasta"),
    log:
        log_path("05_alignment.log"),
    script:
        "lib/StepAlignment.py"


rule tree:
    input:
        data_path("data_alignment.pkl"),
        out_path("clustiso.fasta"),
    output:
        data_path("data_tree.pkl"),
        data_path("dAlTree_tree.pkl"),
        out_path("tree.phylip"),
        out_path("tree.phylip_phyml_tree.txt"),
        out_path("tree.phylip_phyml_stats.txt"),
    log:
        log_path("06_tree.log"),
    script:
        "lib/StepTree.py"


rule duplication:
    input:
        data_path("data_tree.pkl"),
        data_path("dAlTree_tree.pkl"),
        out_path("tree.phylip"),
        out_path("tree.phylip_phyml_tree.txt"),
        out_path("tree.phylip_phyml_stats.txt"),
    output:
        data_path("data_duplication.pkl"),
        data_path("dAlTree_duplication.pkl"),
        out_path("tree.phylip_phyml_tree.txt_recs.nhx"),
        out_path("tree.phylip_phyml_tree.txt_recs.svg"),
    log:
        log_path("07_duplication.log"),
    script:
        "lib/StepDuplication.py"
