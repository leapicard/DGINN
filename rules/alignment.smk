# rule for nucleotide alignment by MAFFT : fasta -> fasta and check coverage
rule align_mafft:
    input:
        resDir+"align_input_{filename}.fasta"
    output:
        resDir+"align_{filename}_mafft.fasta"
    params:
        config = config_file,
        out_mafft = resDir+"align_{filename}_ORFs_al_mafft.fasta"
    message: "\nStarting mafft nucleotide alignment. Aligning sequences from file {input}, writing alignment in file {output}"
    run:
        shell("mafft --auto --quiet {input} > {params.out_mafft}"),
        shell("python3 scripts/covAln.py {params.out_mafft} {output} {params.config}")


# rule for codon alignment by macse or by prank : fasta -> fasta
rule first_align_codon:
    input:
        resDir+"align_{filename}_mafft.fasta" if parameters["align_nt"] else resDir+"align_input_{filename}.fasta"
    output:
        resDir+"align_{filename}_clustiso.fasta"
    params:
        config = config_file,
        extension = parameters["align_codon"],
        out_codon_param = resDir+"align_{filename}_alcodon",
        macse_param = "-prog refineAlignment -align" if parameters["align_nt"] else "-prog alignSequences -seq"
    message: "\nStarting first codon alignment with {params.extension} from file {input}, writing results in file {output}" \
                "\nStarting codon alignment with {params.extension} from file {input}, writing results in file {output}"
    run:
        if parameters["align_codon"][0] == "prank" :
            shell("prank -d={input} -o={params.out_codon_param}_{params.extension}.fas -codon -F"),
            shell("python3 scripts/clusterIso.py {params.config} {output} {params.out_codon_param}_{params.extension}.fas.best.fas")

        elif parameters["align_codon"][0] == "macse" :
            shell("java -jar macse_v2.06.jar {params.macse_param} {input} -out_NT {params.out_codon_param}_{params.extension}.fas"),
            shell("python3 scripts/clusterIso.py {params.config} {output} {params.out_codon_param}_{params.extension}.fas")
        elif len(parameter["align_codon"])==0 and parameters["align_nt"] :
            shell("mv results/align_{filename}_mafft.fasta results/align_{filename}_clustio.fasta")
        else :
            print(f"The codon_aligner parameter in the config file has not been properly filled. The only two values it can contain are 'macse' or 'prank'.It currently contains the value {parameters['align_codon']}")


rule second_align_codon:
    input:
        resDir+"align_{filename}_clustiso.fasta"
    output:
        resDir+"tree_input_{filename}.fasta"
    params:
        config = config_file,
        extension = parameters["align_codon"],
        out_codon_param = resDir+"align_{filename}_alcodon",
        macse_param = "-prog refineAlignment -align" if parameters["align_nt"] else "-prog alignSequences -seq"
    message: "\nStarting second codon alignment with {params.extension} from file {input}, writing results in file {output}"
    run:
        if parameters["align_codon"][0] == "prank":
            shell("mv results/align_{filename}_clustiso.fasta results/tree_input_{filename}.fasta") 
            shell("python3 scripts/clusterIso.py {params.config} {output} {params.out_codon_param}_{params.extension}.fas.best.fas")
        elif parameters["align_codon"][0] == "macse":
            shell("mv results/align_{filename}_clustiso.fasta results/tree_input_{filename}.fasta") 
            shell("python3 scripts/clusterIso.py {params.config} {output} {params.out_codon_param}_{params.extension}.fas")
        elif len(parameter["align_codon"])==0 and parameters["align_nt"] :
            shell("mv results/align_{filename}_clustiso.fasta results/tree_input_{filename}.fasta")
        else :
            print(f"The codon_aligner parameter in the config file has not been properly filled. The only two values it can contain are 'macse' or 'prank'.It currently contains the value {parameters['align_codon']}")
