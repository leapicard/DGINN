
# rule recombination from the script recombination.py : fast, tree -> fasta
rule recombination:
    input:
        resDir + "tree_input_{filename}.fasta" if parameters["step"] == "blast" else resDir + "aln_{filename}.fasta",
        resDir + "tree_{filename}_filtered2.phylip_phyml_tree.txt" if parameters["step"] == "blast" else resDir+treename+".tree"
    output:
        resDir+"recomb_{filename}_files_list.txt"
    params:
        config = config_file
    threads: 6
    message: "\nStarting Tree building, writing output in file {output}" 
    shell:
        "python3 scripts/recombination.py {params.config}"