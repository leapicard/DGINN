# rule blast : enter a fasta file and execute blast from ncbi : fasta -> tsv
rule blast:
    input:
        resDir+"blast_input_{filename}.fasta"
    output:
        resDir+"accessions_input_{filename}.tsv"
    params:
        config = config_file
    message: "\nStarting BLAST alignment on file {input}, writing output in file {output}"
    shell:
        "python3 scripts/blast.py {params.config}"

# rule extract which need NCBI tabulated format :  tsv -> txt
#           and from a list of accession identifiers : txt -> fasta
rule fastaExtract:
    input:
        resDir+"accessions_input_{filename}.tsv" 
    output:
        extract = resDir+"fasta_input_{filename}.txt" ,
        fastaResult = resDir+"orf_input_{filename}.fasta"
    params:
        config = config_file
    message: "\nStarting Accessions search from file {input}, writing accessions in file {output.extract} and {output.fastaResult}"
    shell:
        "python3 scripts/extract.py {params.config} {output.extract} {output.fastaResult}"
        

# rule get orf form mRNA sequences of orthologs : fasta -> fasta
rule getOrf:
    input:
        resDir+"orf_input_{filename}.fasta" 
    output:
        allORFs = resDir+"orf_{filename}_allORFs.fasta",
        orf_longest = resDir+"align_input_{filename}.fasta"
    params:
        config = config_file
    message: "\nStarting ORFs search in FASTA sequences from file {input}.\nWriting longest ORFs found in file {output}"
    shell:
        "python3 scripts/orf.py {params.config} {output.allORFs} {output.orf_longest}"

