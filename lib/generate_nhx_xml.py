import json, sys, re
from Bio import Phylo
from genere_speciestree_xml import createPhyloXML

def phylotree(path_tree):
    treefile = open(path_tree,"r")
    
    folder = path_tree.split(".")[0]
    result = ""
    if path_tree.endswith("phylip_phyml_tree.txt"):
        for line in treefile:
            print(line)
            tline = re.split(' ',line)
            newick=tline[0]
            result = result + createPhyloXML(newick) + "\n"
    elif path_tree.endswith(".nhx"):
        Phylo.convert(path_tree, "newick", folder+".xml", "phyloxml")
        result = Phylo.read(folder+".xml","phyloxml")
    return result

path_tree = "/home/etudiant/Documents/DGINN-API/dginn-api-stage/DGINN/results/ex_aln_recs2.nhx"
print(phylotree(path_tree))

