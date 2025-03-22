import sys
import logging, os
import ete3
from ete3 import PhyloTree
import subprocess, shutil, random
import FastaResFunc, AnalysisFunc
from Bio import SeqIO, Phylo
from statistics import median, mean
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator, DistanceMatrix
from Bio.Phylo import Newick

"""
File which countain all functions about treerecs and tree treatement.
"""

def splitTree(parameters, step="duplication"):
  """
  Split the gene tree in several sub-trees, following several methods.

  1- Cutlongbranches
  2- Reconciliation

  @output The dictionnary {new query, new subalignment file}
  """

  nbspecies=parameters["nbspecies"]

  aln = parameters["input"].split()[0].strip()
  tree = parameters["input"].split()[1].strip()
  outdir = parameters["outdir"]
  
  logger = logging.getLogger(".".join(["main", step]))

  dSubAln = cutLongBranches(parameters, aln, tree, nbspecies, logger)

  dqaln={}

  if parameters["sptree"]!="":
    sptree = parameters["sptree"]
    logger.info("Species tree " + sptree)
    for query, aln in dSubAln.items():
      logger.info("Running Treerecs for " + query)
      recTree = runTreerecs(query, aln, tree, sptree, outdir, logger)
      if recTree:
        dqaln.update(treeParsing(query, aln, recTree, nbspecies, outdir, logger))
      else:
        dqaln[query]=aln
          
  else:
    dqaln=dSubAln

  return(dqaln)


def cutLongBranches(parameters, aln, tree, nbSp, logger):
    """
    Check for overly long branches in a tree and separate both tree and corresponding alignment if found.

    @param1 aln: Fasta alignment
    @param2 tree: Tree corresponding to the alignment
    @param3 logger: Logging object
    @return dSubAln: Updated dictionary of queries and their corresponding alignment file
    """

    LBOpt = parameters["LBopt"]
    
    logger.info("Looking for long branches.")
    loadTree = ete3.Tree(tree)
    dist = [leaf.dist for leaf in loadTree.traverse()]
    # longDist = 500
    dSubAln={}

    queryName=parameters["queryName"]
    outdir=parameters["outdir"]
    
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
            gp = [node] + node.get_children()
            lNewGp = node.get_leaf_names()

            dNewAln = {gene: dID2Seq[gene] for gene in lNewGp if gene in dID2Seq}

            for k in dNewAln:
                dID2Seq.pop(k, None)

            # create new file of sequences

            if len(dNewAln) > nbSp - 1:
              newQuery = queryName +  "_part" + str(matches.index(node) + 1)
              alnf = outdir + "/" + newQuery + ".fasta"
              with open(alnf,"w") as fasta:
                fasta.write(FastaResFunc.dict2fasta(dNewAln))
                fasta.close()
              dSubAln[newQuery] = alnf
            elif len(dNewAln)!=0 or len(dID2Seq)==0:
                logger.info(
                    "Sequences {} will not be considered for downstream analyses as they do not compose a large enough group.".format(
                        " ".join(dNewAln.keys())
                    )
                )

        newQuery = queryName + "_part" + str(len(matches) + 1) 
        alnLeft = os.path.join(outdir,newQuery + ".fasta")

        if len(dID2Seq) > nbSp - 1:
            with open(alnLeft, "w") as fasta:
                fasta.write(FastaResFunc.dict2fasta(dID2Seq))
                logger.info("\tNew alignment:%s" % {alnLeft})
                fasta.close()
            dSubAln[newQuery]=alnLeft
        elif len(dID2Seq)!=0 or len(dID2Seq)==0:
            logger.info(
                "Sequences in {} will not be considered for downstream analyses as they do not compose a large enough group.".format(
                  " ".join(dID2Seq.keys())
                )
            )

    else:
      logger.info("No long branches found.")
      dSubAln[queryName]=aln
      
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


def buildSpeciesTree(sptree, gfaln):
    """
    Build a species tree from the ncbi taxonomy or a given species tree,
    and write it in a specific species tree file.

    @param1 sptree: general species tree.
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
    
    spTreeFile = "/".join(gfaln.split("/")[:-1] + ["species_tree.tree"])
    sptree.write(format=9, outfile=spTreeFile)

    return spTreeFile


def treeCheck(treePath, alnf, logger):
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
      treePath = buildSpeciesTree(treePath, alnf)

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
    
    @return dquery: dictionnary (new queries, new files)
    """

    with open(recTree, "r") as tree:
        reconTree = tree.readlines()[1]
        tree.close()

    testTree = ete3.Tree(reconTree)

    seqs = SeqIO.parse(open(ORF), "fasta")
    dID2Seq = {gene.id: gene.seq for gene in seqs}

    # get all nodes annotated with a duplication event
    dupl = testTree.search_nodes(D="Y")
    dNb2Node = {int(node.ND): node for node in dupl}
    nDuplSign = 0
    dOut = {}
    sp = set([leaf.S for leaf in testTree])
    dDupl2Seq = {}

    # as long as the number of species left in the tree is equal or superior to the cut-off specified by the user and there still are nodes annoted with duplication events
    while len(sp) > int(nbSp) - 1 and len(dNb2Node.keys()) > 0:
        # start from the most recent duplications (ie, the furthest node)
        sp = set([leaf.S for leaf in testTree])
        nodeNb = min(dNb2Node.keys())
        node = dNb2Node[nodeNb]

        # for each of the branches concerned by the duplication
        nGp = 1
        interok = False
        
        # do not consider dubious duplications (no intersection between the species on either side of the annotated duplication)
        lf = [set([leaf.S for leaf in gp]) for gp in node.get_children()]
        interok = (
            len(lf[0].intersection(lf[1])) != 0
            and len(lf[0]) > int(nbSp) / 2 - 1
            and len(lf[1]) > int(nbSp) / 2 - 1
        )

        if not interok:
            dNb2Node.pop(nodeNb, None)

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
                        newQuery = query + "_D%d_gp%d"%(nodeNb,nGp)
                        outFile = os.path.join(outdir, newQuery + "_orf.fasta")
                        dOut[newQuery]=outFile

                        # create new file of orthologous sequences
                        with open(outFile, "w") as fasta:
                            fasta.write(FastaResFunc.dict2fasta(dOrtho2Seq))
                            fasta.close()
                        # remove the node from the tree
                        removed = gp.detach()

                    dDupl2Seq["{:d}-{:d}".format(nodeNb, nGp)] = orthos
                nGp += 1

            dNb2Node.pop(nodeNb, None)

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
        dOut[newQuery]=outFile
      else:
        logger.info(
          "Ignoring remaining sequences {} as they do not compose a group of enough orthologs.".format(
            list(dRemain.keys())
          )
        )

    # check that all files contain sequences, otherwise filter them out
    for dupKey, dupFile in dOut.items():
        lnseq = len([seq for seq in SeqIO.parse(open(dupFile), "fasta")])
        if lnseq < nbSp:
          dOut.remove(dupKey)
          os.remove(dupFile)

    logger.info(
        "{:d} duplications detected by Treerecs, extracting {:d} groups of at least {} orthologs.".format(
            len(dupl), len(dOut), nbSp
        )
    )
    return dOut


###
# Post order traversal of sptree to find polytomies


def max_parcimony_polytomy(node, gtree, outdir):
    ch = node.get_children()
    lch = len(node.get_children())
    # best newtree
    bnt = False

    lspt = []
    lcost = []

    otmp=os.path.join(outdir,"tmp")
    if not os.path.exists(otmp):
        os.makedirs(otmp)

    gf = os.path.join(otmp, "eval_poly_gene.tree")
    gtree.write(outfile=gf)
    
    for desctree in enum_unordered(range(lch)):
        nt = PhyloTree(str(desctree) + ";")
        for p in nt.get_leaves():
            if int(p.name) == lch:
                continue

            p.add_sister(PhyloTree(ch[int(p.name)].write()))
            p.detach()

        lspt.append(nt)
        spf = os.path.join(otmp, "eval_poly_sp_%d.tree" % (len(lspt)))
        nt.write(outfile=spf)
        gfi = os.path.join(otmp, "eval_poly_gene_%d.tree" % (len(lspt)))
        if os.path.exists(gfi):
          os.remove(gfi)

        os.symlink(gf, gfi)

        val = "treerecs -g {:s} -s {:s} -o {:s} -f -t 0.8 -O NHX".format(gfi, spf, otmp)
        subprocess.run(val, shell=True, capture_output=True)

        fnt = open(os.path.join(otmp, "eval_poly_gene_%d.tree_recs.nhx" % (len(lspt))), "r")
        lc = fnt.readline()
        fnt.close()
        pc = lc.find("total cost")
        pe = lc.find("=", pc)
        lcost.append(int(lc[pe + 1 : lc.find(",", pc)]))

    # now keep the most parcimonious

    m = min(lcost)
    im = [i for i in range(len(lcost)) if lcost[i] == m]

    # clean tmp
    shutil.rmtree(otmp)

    return lspt[im[0]]  # get first reconciliation if equality...


###
# Resolve high order polytomy through nj on sampled clades

def nj_sample_polytomy(node, gtree, outdir):
    """ Resolve polytomy from observed gene tree."""

    def g2sp(gname):
      return("_".join(gname.split("_")[:2]))
      
    ## get list of lists of leaves names for all children
    lln=[]
    lch = node.get_children()
    for ch in lch:
      lln.append(ch.get_leaf_names())

    ## set genes leaves names as species names
    obstree=PhyloTree(gtree.write())
    
    ### Sample quartets of species names to perform partial reconciliations
    nbSample = 50
    nbleaves = 5
    
    distmat=[]
    cummat=[]
    for i in range(len(lch)):
      distmat.append([0]*(i+1))
      cummat.append([0]*(i+1))
   
    lnoderes=[]
    for ns in range(nbSample):
      inodes = random.sample(range(len(lch)), nbleaves)
      snodes = [lch[i] for i in inodes]
      
      node2 = PhyloTree("Root;")
      for ch in snodes:
        ch2= PhyloTree(ch.write())
        node2.add_child(ch2)

      gt2 = gtree.copy()
      leavnode = [l.name for l in node2.get_leaves()]
      leag = [l.name for l in gt2.get_leaves() if g2sp(l.name) in leavnode]
      if len(leag) != 0:
        gt2.prune(leag)

      noderes=max_parcimony_polytomy(node2, gt2, outdir)
      lnoderes.append(noderes)

      for isa in inodes:
        for jsa in inodes:
          if isa<=jsa:
            continue
          d=noderes.get_distance(lch[isa].get_leaf_names()[0],lch[jsa].get_leaf_names()[0],topology_only=True)
          distmat[isa][jsa]+=d
          cummat[isa][jsa]+=1
      
    ### Compute mean distances between all leaves from samples

    def rap(x,y):
      if y==0:
        return 0
      else:
        return x/y
      
    meanmat = [[rap(distmat[i][j],cummat[i][j]) for j in range(len(distmat[i]))] for i in range(len(distmat))]
    
    dm=DistanceMatrix(names=list(map(str,range(len(lch)))),matrix=meanmat)

    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    tmpt = os.path.join(outdir,"tmp.dnd")
    Phylo.write(tree,tmpt,format="newick") # do not know how to make otherwise than through a temp file

    with open(tmpt, "r") as tree:
        restree= tree.readlines()[0]
        tree.close()
    # clean tmp file
    os.remove(tmpt)
    treeok = PhyloTree(restree,format=1)

    ## Put back correct children
    for num in range(len(lch)):
      chn = treeok&str(num)
      chn.add_sister(lch[num])
      chn.detach()
      
    return treeok


def runTreerecs(query, aln, pathGtree, pathSptree, outdir,logger):
    """
    Procedure which launch treerecs. 

    @param1 query: name of the alignment
    @param2 aln: file name of the alignment
    @param1 pathGtree: Path to the gene tree file
    @param2 pathSptree: Path to the species tree
    @param3 outdir: Output directory

    @output file name of reconciliated tree
    """

    ## prune gene tree according to species Tree
    try:
      pathGtree = filterTree(pathGtree, pathSptree)
    except TreeError:
      ## names of the genes
      seqs = SeqIO.parse(open(aln), "fasta")
      dID2Seq = [gene.id for gene in seqs]
      ## prune gene tree according to aln sequences
      pathGtree = PhyloTree(pathGtree)
      pathGtree.prune(dID2Seq)

    ## look for polytomies, and change species tree in a most
    ## parcimonious way

    gtree = PhyloTree(pathGtree)
    sptree = PhyloTree(pathSptree)
  
    # species of the genes
    lg = gtree.get_leaf_names()
    gs = set(["_".join(g.split("_")[:2]) for g in lg])

    ## prune species tree
    sptree.prune(gs)
    
    ## look for polytomies, and change species tree in a most
    ## parcimonious way for species under the polytomy (but not with
    ## the others!).

    thrspoly = 6
    poly = False

    lnode = [node for node in sptree.traverse("postorder")]
    for node in lnode:
        lch = len(node.get_children())
        while lch >2:
          if not poly:
            logger = logging.getLogger("main.duplication")
            logger.warning(
              "Species Tree with polytomies: solved with most parcimonious or neighbour-joining."
            )
            poly = True

          gt2 = gtree.copy()
          leavnode = [l.name for l in node.get_leaves()]
          leag = [l.name for l in gt2.get_leaves() if "_".join(l.name.split("_")[:2]) in leavnode]
          if len(leag) != 0:
            gt2.prune(leag)

          if lch> thrspoly:
            nb = nj_sample_polytomy(node, gt2, outdir)
          else:
            nb = max_parcimony_polytomy(node, gt2, outdir)

          node.add_sister(nb)
          node.detach()
          node=nb
          lch = len(node.get_children())

    if not poly:  # no polytomy solved
        pathSptree2 = pathSptree
    else:
        lp = pathSptree.split(".")

        pathSptree2 = ".".join(lp[:-1] + ["bif"] + [lp[-1]])
        sptree.write(outfile=pathSptree2)

    ### filter out unmatched genes in species tree
    val = "treerecs -g {:s} -s {:s} -o {:s} -f -t 0.8 -O NHX:svg".format(
        pathGtree, pathSptree2, outdir
    )

    subprocess.run(val, shell=True, stdout = subprocess.PIPE)
    
    return os.path.join(pathGtree + "_recs.nhx")  # o+suff[:suff.rfind(".")]+"_recs.nhx"




#######=================================================================================================================
