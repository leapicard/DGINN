
# rule the first phylogenetic tree : fasta -> txt

rule tree:
    input:
        resDir+"tree_input_{filename}.fasta"
    output:
        resDir+"tree_{filename}_filtered2.phylip_phyml_tree.txt"
    params:
        fonction = "phyMLTree",
        config = config_file
    message: "\nStarting Tree building, writing output in file {output}" 
    shell:
        "python3 scripts/tree.py {params.config} {params.fonction} {output}"
