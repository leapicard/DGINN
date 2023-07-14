

# rule for calculate positive selection from positiveSelection.py

rule positiveSelection:
    input:
        resDir+"tree_{filename}_filtered2.phylip_phyml_tree.txt" if (parameters["step"] == "blast" and not parameters["duplication"])
        else resDir+"dup_{filename}_recs2.nhx" if (parameters["step"] == "blast" and parameters["duplication"]) or (parameters["step"] == "option" and parameters["duplication"])
        else resDir+treename+".tree"
    output:
        resDir+"possel_{filename}_files_list.txt" 
    threads: 8
    params:
        config = config_file
    run:
        shell("python3 scripts/positiveSelection.py {params.config}")