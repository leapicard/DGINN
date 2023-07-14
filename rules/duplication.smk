


# rule duplication from the script Duplication.py : fast, tree -> nhx
rule duplication:
    input:
        resDir+"tree_input_{filename}.fasta" if parameters["step"]=="blast" else resDir+"aln_{filename}.fasta",
        resDir+"tree_{filename}_filtered2.phylip_phyml_tree.txt" if parameters["step"]=="blast" else resDir+treename+".tree"
    params:
        config = config_file
    output:
        resDir+"dup_{filename}_recs2.nhx"
    threads: 6
    message: "\nStarting Tree building, writing output in file {output}" 
    shell:
        "python3 scripts/duplication.py {params.config} {output}"