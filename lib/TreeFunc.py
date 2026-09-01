import sys
import logging, os, re
import ete3
from ete3 import PhyloTree
import subprocess, shutil, random
import FastaResFunc, AnalysisFunc
from Bio import SeqIO, Phylo
from statistics import median, mean
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator, DistanceMatrix
from Bio.Phylo import Newick
from ete3.coretype.tree import TreeError

"""
File which countain all functions about treerecs and tree treatement.
"""

def splitTree(parameters, step="duplication"):
  """Split the gene tree in several sub-trees, following several methods.

  1- Reconciliation (with re-rooting)
  2- Cutlongbranches
  3- Reduce polymorphism
  
  @output The dictionnary {new query, new subalignment file}. If
  nothing new is built returns {query, alignment file}

  """

  nbspecies=parameters["nbspecies"]
  
  aln = parameters["input"].split()[0].strip()
  tree = parameters["input"].split()[1].strip()
  outdir = parameters["outdir"]
  query = parameters["queryName"]
  poly = parameters["SNP"]
  
  logger = logging.getLogger(".".join(["main", step]))

  ### Reconciliation
  dqaln={}
  if parameters["sptree"]!="":
    sptree = parameters["sptree"]
    logger.info("Species tree " + sptree)
    recTree = runNotung(query, aln, tree, sptree, outdir, logger)
    if recTree:
      dqaln.update(treeParsing(query, aln, recTree, nbspecies, outdir, logger))
    else:
      dqaln[query]=[aln, tree]

  if len(dqaln)==0:
    dqaln[query]=[aln, tree]

    
  ### cutLongBranches

  dSubAln = {}
  
  for query, [aln, tree] in dqaln.items():
    dSubAln.update(cutLongBranches(query, aln, tree, parameters, nbspecies, outdir, logger))

  ### merge polymorphism
  if poly:
    dqaln = {}
    for query, [aln, tree] in dSubAln.items():
      logger.info("merge "+ query)
      k, aln = AnalysisFunc.mergePolymorphism(query, aln, tree, outdir, logger)
      dqaln[k] = aln
  else:
    dqaln = {q:aln for q,[aln,tree] in dSubAln.items()}
      
  return(dqaln)


def cutLongBranches(queryName, aln, tree, parameters, nbSp, outdir, logger):
    """
    Check for overly long branches in a tree and separate both tree and corresponding alignment if found.

    @param1 queryName: query
    @param1 aln: Fasta alignment
    @param2 tree: Tree corresponding to the alignment
    @param3 parameters: used parameters for cutLongBranches
    @param4 outdir: output directory
    @param3 logger: Logging object
    @return Dictionary of {queries,[alignment file, tree file]}
    """

    LBOpt = parameters["LBopt"]
    
    logger.info("Looking for long branches.")

    loadTree = ete3.Tree(tree)
    dist = [leaf.dist for leaf in loadTree.traverse()]
    # longDist = 500
    dSubAln={}

    if "cutoff" in LBOpt:
        if "(" in LBOpt:
            factor = float(LBOpt.split("(")[1].replace(")", ""))
        else:
            factor = 50
        medianDist = median(dist)
        meanDist = mean(dist)
        longDist = meanDist * factor
    elif "IQR" in LBOpt:
        if "(" in LBOpt:
            factor = float(LBOpt.split("(")[1].replace(")", ""))
        else:
            factor = 50
        df = pd.DataFrame(dist)
        Q1 = df.quantile(0.25)
        Q3 = df.quantile(0.75)
        IQR = Q3 - Q1
        lDist = Q3 + (factor * IQR)
        longDist = lDist[0]

    logger.info(
        "Long branches will be evaluated through the {} method (factor {})".format(
            LBOpt, factor
        )
    )
    nbSp = int(nbSp)
    matches = [leaf for leaf in loadTree.traverse("postorder") if leaf.dist > longDist]
    if len(matches) > 0:
        logger.info(
            "{} long branches found, separating alignments.".format(len(matches))
        )

        seqs = SeqIO.parse(open(aln), "fasta")
        dID2Seq = {gene.id: gene.seq for gene in seqs}

        for node in matches:
            up = node.up
            gp = node.detach()
            lNewGp = gp.get_leaf_names()

            # iteratively remove nodes with on child
            while up and len(up.children)==1:
              upup = up.up
              up.delete()
              up = upup
              
            dNewAln = {gene: dID2Seq[gene] for gene in lNewGp if gene in dID2Seq}

            for k in dNewAln:
                dID2Seq.pop(k, None)

            # create new file of sequences

            if len(dNewAln) > nbSp - 1:
              newQuery = queryName +  "_part" + str(matches.index(node) + 1)
              alnf = outdir + "/" + newQuery + "_orf.fasta"
              with open(alnf,"w") as fasta:
                fasta.write(FastaResFunc.dict2fasta(dNewAln))
                fasta.close()
              outTree = outdir + "/" + newQuery + "_orf.dnd"
              gp.write(outfile = outTree)
              dSubAln[newQuery] = [alnf,outTree]
            elif len(dNewAln)!=0 or len(dID2Seq)==0:
                logger.info(
                    "Sequences {} will not be considered for downstream analyses as they do not compose a large enough group.".format(
                        " ".join(dNewAln.keys())
                    )
                )

        newQuery = queryName + "_part" + str(len(matches) + 1) 
        alnLeft = os.path.join(outdir,newQuery + "_orf.fasta")
        treeLeft = os.path.join(outdir,newQuery + "_orf.dnd")

        if len(dID2Seq) > nbSp - 1:
            with open(alnLeft, "w") as fasta:
                fasta.write(FastaResFunc.dict2fasta(dID2Seq))
                logger.info("\tNew alignment:%s" % {alnLeft})
                fasta.close()
                
            loadTree.write(outfile=treeLeft)
            dSubAln[newQuery]=[alnLeft,treeLeft]


        elif len(dID2Seq)!=0 or len(dID2Seq)==0:
            logger.info(
                "Sequences in {} will not be considered for downstream analyses as they do not compose a large enough group.".format(
                  " ".join(dID2Seq.keys())
                )
            )

    else:
      logger.info("No long branches found.")
      dSubAln[queryName]=[aln,tree]
      
    return dSubAln

#######################################
#### Class used for resolving polytomies (from Stackoverflow)


# A very simple representation for Nodes. Leaves are anything which is not a Node.
class Node(object):
    def __init__(self, left, right):
        self.left = left
        self.right = right

    def __repr__(self):
        return "(%s, %s)" % (self.left, self.right)


# Given a tree and a label, yields every possible augmentation of the tree by
# adding a new node with the label as a child "above" some existing Node or Leaf.
def add_leaf(tree, label):
    yield Node(label, tree)
    if isinstance(tree, Node):
        for left in add_leaf(tree.left, label):
            yield Node(left, tree.right)
        for right in add_leaf(tree.right, label):
            yield Node(tree.left, right)


# Given a list of labels, yield each rooted, unordered full binary tree with
# the specified labels.
def enum_unordered(labels):
    if len(labels) == 1:
        yield labels[0]
    else:
        for tree in enum_unordered(labels[1:]):
            for new_tree in add_leaf(tree, labels[0]):
                yield new_tree


#######Treerecs=========================================================================================================
# =========================================================================================================================
def getLeaves(path):
    """
    Open a newick file and return a list of species

    @param path: path of a tree file
    @return lTreeData: list of species
    """
    tree = ete3.Tree(path)
    lGene = tree.get_leaf_names()

    return lGene
 

def buildSpeciesTree(queryName, gfaln):
    """
    Build a species tree from the ncbi taxonomy or a given species tree,
    and write it in a specific species tree file.

    @param1 queryName: query name.
    @param2 gfaln: path of the alignment
    
    @return the species tree file path.
    """

    
    ncbi = ete3.NCBITaxa(dbfile="/opt/ncbi/taxa.sqlite")
    
    ## get names from link Taxids & species names abbreviations

    accns = list(SeqIO.parse(gfaln, "fasta"))
    gLeavesSp = [" ".join(name.id.split("_")[:2]) for name in accns]

    gLeavesid = ncbi.get_name_translator(gLeavesSp)
    sptree = ncbi.get_topology([v[0] for v in gLeavesid.values()])

    for node in sptree.traverse():
      node.name = ncbi.get_taxid_translator([int(node.name)])[int(node.name)]
      node.name = "_".join(node.name.split())

    spTreeFile = "/".join(gfaln.split("/")[:-1] + [queryName+"_species_tree.tree"])
    sptree.write(format=9, outfile=spTreeFile)

    return spTreeFile


def treeCheck(treePath, alnf, queryName, logger):
    """
    Check if the tree isn't corrupted

    @param1 treePath: tree's path
    @param2 alnf: alignment's path
    @return1 treePath: tree's path
    """

    buildspt = False
    
    if treePath != "":
      if not os.path.exists(treePath):
        logger.warning("Species tree file " + treePath + " does not exist.")
        buildspt = True
      elif not ete3.Tree(treePath):
        logger.warning("The species tree is corrupted.")
        buildspt = True
    else:
      buildspt = True

    if buildspt:
      treePath = buildSpeciesTree(queryName, alnf)

    return treePath


# =========================================================================================================================


def assocFile(sptree, path, dirName):
    """
    Create a file which contain the species of each genes

    @param1 sptree: Species tree
    @param2 path: Path of a fasta file
    @param3 dirName: Name for a new directory
    """

    lGeneId = []
    ff = open(path, "r")
    for accn in SeqIO.parse(ff, "fasta"):
        lGeneId.append(accn.id)
    ff.close()
    lGeneId.sort()

    lTreeData = getLeaves(sptree)
    lTreeData.sort()
    dSp2Gen = {}
    index = 0
    # we go through the list one and two in the same time to gain time to associate genes with their species.
    for sp in lTreeData:
        indexSave = index
        while sp.lower() not in lGeneId[index].lower():
            if index + 1 < len(lGeneId):
                index += 1
            else:
                index = indexSave
                break

        while sp.lower() in lGeneId[index].lower():
            dSp2Gen[lGeneId[index]] = sp
            if index + 1 < len(lGeneId):
                index += 1
            else:
                break

    out = dirName + "/" + path.split("/")[-1].split(".")[0] + "_Species2Genes.txt"

    with open(out, "w") as cor:
        for k, v in dSp2Gen.items():
            cor.write(k + "\t" + v + "\n")
        cor.close()

    return out


def supData(filePath, corFile, dirName):
    """
    Delete genes in a fasta file

    @param1 filePath: Path to a fasta file
    @param2 corFile: Path to the file with correspondence beetween gene and species
    @param3 dirName: Name of a directory
    @return out: Path
    """

    with open(corFile, "r") as corSG:
        lcorSG = corSG.readlines()
        lGene = [i.split("\t")[0] for i in lcorSG]
        corSG.close()
    newDico = {}
    for accn in SeqIO.parse(open(filePath, "r"), "fasta"):
        if accn.id in lGene:
            newDico[accn.id] = accn.seq

    out = dirName + "/" + filePath.replace(".fasta", "_filtered.fasta").split("/")[-1]
    with open(out, "w") as newVer:
        newVer.write(FastaResFunc.dict2fasta(newDico))
        newVer.close()
    return out


# =========================================================================================================================


def filterTree(tree, spTree):
    """
    Delete genes in the tree which aren't in the species tree

    @param1 tree: Path to a tree file
    @param2 spTree: Path to a species tree file
    @return out: Path to the filtered tree
    """

    tg = ete3.Tree(tree)

    ts = ete3.Tree(spTree)
    lg = tg.get_leaf_names()
    ls = ts.get_leaf_names()

    #gok = [g for g in lg if g.split("_")[] in ls]
    gok = [g for g in lg if "_".join(g.split("_")[:2]) in ls]

    tg.prune(gok)
    out = tree.replace(".tree", "_filtered.tree")
    tg.write(outfile=out)

    return out


# =========================================================================================================================
def treeParsing(query, ORF, recTree, nbSp, outdir, logger):
    """
    Function which parse gene data in many group according to duplication in the reconciliated tree

    @param1 query 
    @param2 ORFs: Path to the ORFs file
    @param3 recTree: Path to a reconciliation tree file
    @param4 nbSp: threshold number of leaves to build a separate clade
    @param5 outdir: Output directory
    @param6 logger: Logging object
    
    @return dquery: dictionnary {new queries, [new alignment, subtree]}
    """

     
    with open(recTree, "r") as tree:
      reconTree = tree.readlines()[0]
      tree.close()

    testTree = ete3.PhyloTree(reconTree, format=1)

    seqs = SeqIO.parse(open(ORF), "fasta")
    dID2Seq = {gene.id: gene.seq for gene in seqs}
    
    # get all nodes annotated with a duplication event, in pre-order strategy
    
    dupl = testTree.search_nodes(D="Y")
    dNb2Node = [node for node in dupl]
    dNb2Node.reverse()

    nDuplSign = 0
    dOut = {}
    sp = set([leaf.S for leaf in testTree])
    dDupl2Seq = {}

    # as long as the number of species left in the tree is equal or superior to the cut-off specified by the user and there still are nodes annoted with duplication events
    while len(sp) > int(nbSp) - 1 and len(dNb2Node) > 0:
        # start from the most recent duplications (ie, the furthest node)
        sp = set([leaf.S for leaf in testTree])
        node = dNb2Node[0]

        # for each of the branches concerned by the duplication
        nGp = 1

        ###################
        ### here deactivated
        # do not consider dubious duplications (no intersection between the species on either side of the annotated duplication)
        interok = False
        lf = [set([leaf.S for leaf in gp]) for gp in node.get_children()]
        interok = (
            len(lf) > 1 and 
            len(lf[0].intersection(lf[1])) != 0
            and len(lf[0]) > int(nbSp) / 2 - 1
            and len(lf[1]) > int(nbSp) / 2 - 1
        )

        #if not interok: 
        #    dNb2Node.pop(0)

        if False:
          pass
        ####################
        # otherwise check it out
        else:
            for gp in node.get_children():
                spGp = set([leaf.S for leaf in gp])

                # check if the numbers of species in the branch is equal or superior to the cut-off specified by the user
                if len(spGp) > int(nbSp) - 1:
                    orthos = gp.get_leaf_names()
                    dOrtho2Seq = {
                        ortho: dID2Seq[ortho]
                        for ortho in orthos
                        if not ortho == "" and ortho in dID2Seq
                    }

                    # check if orthologs have already been included in another, more recent, duplication event
                    already = False
                    for doneDupl in dDupl2Seq:
                        if all(ortho in dDupl2Seq[doneDupl] for ortho in orthos):
                            already = True
                            break

                    if not already:
                        nDuplSign += 1
                        newQuery = query + "_D%s_gp%d"%(node.name[1:],nGp)
                        outFile = os.path.join(outdir, newQuery + "_orf.fasta")

                        # create new file of orthologous sequences
                        with open(outFile, "w") as fasta:
                            fasta.write(FastaResFunc.dict2fasta(dOrtho2Seq))
                            fasta.close()
                        # remove the node from the tree
                        parent = gp.up
                        removed = gp.detach()

                        # clean empty branches
                        while parent is not None:
                          grand_parent = parent.up
                          if grand_parent is None:
                            break
                          if len(parent.children) == 0: # if no child, remove
                            parent.detach()
                          elif len(grand_parent.children) == 1: # if on false node, plug children to grand parent 
                            for ch in parent.children:
                              grand_parent.add_child(ch)
                          else:  ## no need to go up again
                            break
                          parent = grand_parent

                        ## write detached tree
                        outTree = os.path.join(outdir, newQuery + "_orf.dnd")
                        removed.write(outfile=outTree)
                        dOut[newQuery]=[outFile,outTree]
                        
                    logger.info("Extracting clade of {:d} species under node {:s}".format(len(spGp),node.name[1:]))

                    dDupl2Seq["{:s}-{:d}".format(node.name[1:], nGp)] = orthos
                nGp += 1

            dNb2Node.pop(0)

    # if duplication groups have been extracted
    # pool remaining sequences (if span enough different species - per user's specification) into new file
    if len(dOut) > 0:
      leftovers = filter(None, testTree.get_leaf_names())
      dRemain = {left: dID2Seq[left] for left in leftovers if left in dID2Seq}

      if len(dRemain.keys()) > int(nbSp) - 1:
        newQuery = query + "_Drem"
        outFile = os.path.join(outdir, newQuery + "_orf.fasta")
        nDuplSign += 1
        with open(outFile, "w") as fasta:
          fasta.write(FastaResFunc.dict2fasta(dRemain))
          fasta.close()
        outTree = os.path.join(outdir, newQuery + "_orf.dnd")
        testTree.write(outfile=outTree)
        
        dOut[newQuery]=[outFile, outTree]
        logger.info("Extracting remaining sequence of {:d} species".format(len(spGp)))
      else:
        logger.info(
          "Ignoring remaining sequences {} as they do not compose a group of enough orthologs.".format(
            list(dRemain.keys())
          )
        )
                    
    # check that all files contain sequences, otherwise filter them out
    rmKey = []
    for dupKey, [dupFile, dupTree] in dOut.items():
        lnseq = len([seq for seq in SeqIO.parse(open(dupFile), "fasta")])
        if lnseq < nbSp:
          rmKey.append(dupKey)
          os.remove(dupFile)
          os.remove(dupTree)

    for key in rmKey:
      dOut.pop(key)
      
    logger.info(
        "{:d} duplications detected, extracting {:d} groups of at least {} orthologs.".format(
            len(dupl), len(dOut), nbSp
        )
    )
    return dOut



#######=================================================================================================================

def runNotung(query, aln, pathGtree, pathSptree, outdir, logger):
    """
    Procedure which launches Notung. 

    @param1 query: name of the alignment
    @param2 aln: file name of the alignment
    @param1 pathGtree: Path to the gene tree file
    @param2 pathSptree: Path to the species tree
    @param3 outdir: Output directory

    @output file name of reconciliated tree
    """


    ## set arbitrary outgroup at the top of the tree, necessary with species tree with polytomies
    logger.info("run Notung on " + pathGtree)

    gtree = PhyloTree(pathGtree)
    childroot = gtree.get_children()
    if len(childroot)>2:
      gtree.set_outgroup(childroot[0])
    
    gtree.write(format=9, outfile=pathGtree)

    ### pruning, reconciliation & rooting of gene tree

    val = "java -jar lib/Notung-2.9.1.5.jar -s {:s} -g {:s} --prune --root --treeoutput nhx --outputdir {:s} --reconcile --nolosses".format(pathSptree, pathGtree, outdir) # --rearrange --threshold 0.8
    AnalysisFunc.cmd(val,True)

    return os.path.join(outdir,os.path.split(pathGtree)[-1] + ".reconciled")


  ### filter out unmatched genes in species tree
#    val = "treerecs -r -g {:s} -s {:s} -o {:s} -f -t 0.8 -O NHX:svg".format(
#        pathGtree, pathSptree2, outdir
#    )

 
