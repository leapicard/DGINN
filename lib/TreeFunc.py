import logging, os
import ete3
import FastaResFunc, AnalysisFunc, ExtractFunc
from Bio import SeqIO
from collections import defaultdict

"""
File which countain all functions about treerecs and tree treatement.
"""

#######Treerecs=========================================================================================================
#=========================================================================================================================
def getLeaves(path):
	"""
	Open a newick file and return a list of species

	@param path: path of a tree file
	@return lTreeData: list of species
	"""
	tree = ete3.Tree(path)
	lGene = tree.get_leaf_names()

	return lGene


def treeCheck(treePath, optionTree, logger):
	"""
	Check if the tree isn't corruped

	@param1 treePath: tree's path
	@param2 optionTree: Boolean (treerecs option)
	@param3 logger: Logging object
	@return1 treePath: tree's path
	@return2 optionTree: Boolean (treerecs option)
	"""
	if treePath != "":
		if not os.path.exists(treePath):
			logger.warning("The path for the species tree doesn't exist.")
			logger.warning("Treerecs option will be set to False.")
			optionTree = False
		else:
			if not ete3.Tree(treePath):
				logger.warning("The species tree is corrupted.")
				logger.warning("Treerecs option will be set to False.")
				optionTree = False

	return treePath, optionTree

#=========================================================================================================================

def assocFile(sptree, path, dirName):
	"""
	Create a file which contain the species of each genes

	@param1 sptree: Species tree
	@param2 path: Path of a fasta file
	@param3 dirName: Name for a new directory
	"""
		
	lTreeData = getLeaves(sptree)

	lGeneId = []
	for accn in SeqIO.parse(open(path, "r"), "fasta"):
		lGeneId.append(accn.id)

	lGeneId.sort()
	lTreeData.sort()
	dSp2Gen = {}
	index = 0
	#we go through the list one and two in the same time to gain time to associate genes with their species.
	for sp in lTreeData:

		indexSave = index
		while sp.lower() not in lGeneId[index].lower() :
			if index+1 < len(lGeneId):
				index += 1
			else:
				index = indexSave
				break

		while sp.lower() in lGeneId[index].lower() :
			dSp2Gen[lGeneId[index]] = sp
			if index+1 < len(lGeneId):
				index += 1
			else:
				break

	out = dirName+path.split("/")[-1].split(".")[0]+"_Species2Genes.txt"

	with open(out, "w") as cor:
		for i in dSp2Gen:
			cor.write(i+"\t"+dSp2Gen[i]+"\n")

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
		corSG = corSG.readlines()
		lGene = [ i.split("\t")[0] for i in corSG ]

	newDico = {}
	for accn in SeqIO.parse(open(filePath, "r"), "fasta"):
		if accn.id in lGene:
			newDico[accn.id] = accn.seq

	out = dirName+filePath.replace(".fasta", "_filtered.fasta").split("/")[-1]
	with open(out, "w") as newVer:
		newVer.write(FastaResFunc.dict2fasta(newDico))

	return out

#=========================================================================================================================

def filterTree(tree, spTree, cor):
	"""
	Delete genes in the tree which aren't in the species tree

	@param1 tree: Path to a tree file
	@param2 spTree: Path to a species tree file
	@return out: Path
	"""
	t1 = ete3.Tree(tree)
	t2 = ete3.Tree(spTree)
	g1 = t1.get_leaf_names()
	g2 = t2.get_leaf_names()

	lCor = []
	with open(cor, "r") as lCor:
		lCor = lCor.readlines()

	dCor = {}
	dCorInv = {}
	for i in lCor:
		temp = i.split("\t")
		dCor[temp[0]] = temp[1].strip("\n")
		dCorInv[i[1].strip("\n")] = i[0]

	g1Genes = []
	
	for z in g1:
		if z in dCor.keys() and dCor[z] in g2:
			g1Genes.append(z)

	"""if len(g1Genes) < len(g1)/2:
		print(len(g1Genes))
		print(len(g1)/2)
		exit()"""

	t1.prune(set(g1Genes))
	out = tree.replace(".tree", "_filtered.tree")
	t1.write(out)
	
	return out

#=========================================================================================================================
def treeParsing(ORF, recTree, nbSp, o, logger):
	"""
	Function which parse gene data in many group according to duplication in the reconciliated tree

	@param1 ORFs: Path to the ORFs file
	@param2 tree: Path to a tree file
	@param3 geneName: Gene name
	@param4 o: Output directory
	@param5 logger: Logging object
	@return lOut: List of path (fasta files)
	"""
	
	with open(recTree, "r") as tree:
		reconTree = tree.readlines()[1]
		testTree = ete3.Tree(reconTree)
		
		seqs = SeqIO.parse(open(ORF),'fasta')
		dID2Seq = {gene.id: gene.seq for gene in seqs}
		
		# get all nodes annotated with a duplication event
		dupl = testTree.search_nodes(D="Y")
		dNb2Node = {int(node.ND): node for node in dupl}
		nDuplSign = 0
		lOut = []
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
			interok = (len(lf[0].intersection(lf[1])) != 0 and len(lf[0]) > int(nbSp)/2 - 1 and len(lf[1]) > int(nbSp)/2 - 1)
	
			if not interok:
				dNb2Node.pop(nodeNb, None)
			
			# otherwise check it out	
			else:	
				for gp in node.get_children():
					spGp = set([leaf.S for leaf in gp])
					
					# check if the numbers of species in the branch is equal or superior to the cut-off specified by the user
					if len(spGp) > int(nbSp) - 1:
						
						orthos = gp.get_leaf_names()
						dOrtho2Seq = {ortho: dID2Seq[ortho] for ortho in orthos if not ortho == ""}
						
						#check if orthologues have already been included in another, more recent, duplication event
						already = False
						for doneDupl in dDupl2Seq:
							if all(ortho in dDupl2Seq[doneDupl] for ortho in orthos):
								already = True
								break
						
						if not already:
							nDuplSign += 1
							outFile = o+ORF.split("/")[-1].split(".")[0]+"_D"+str(nodeNb)+"gp"+str(nGp)+".fasta"
							lOut.append(outFile)
							
							# create new file of orthologous sequences
							with open(outFile, "w") as fasta:
								fasta.write(FastaResFunc.dict2fasta(dOrtho2Seq))
						
							# remove the node from the tree
							removed = gp.detach()
						
						dDupl2Seq["{:d}-{:d}".format(nodeNb, nGp)] = orthos
					nGp += 1
				
				dNb2Node.pop(nodeNb, None)
			
			# if duplication groups have been extracted
			# pool remaining sequences (if span enough different species - per user's specification) into new file
			if len(lOut) > 0:
				leftovers = filter(None, testTree.get_leaf_names())
				dRemain = {left: dID2Seq[left] for left in leftovers}
				
				if len(dRemain.keys()) > int(nbSp) - 1:
					outFile = o+ORF.split("/")[-1].split(".")[0]+"_duplication_remainingsequences.fasta"
					
					with open(outFile, "w") as fasta:
						fasta.write(FastaResFunc.dict2fasta(dRemain))
					lOut.append(outFile)
				else:
					logger.info("Ignoring remaining sequences {} as they do not compose a group of enough orthologs.".format(list(dRemain.keys())))
				
	logger.info("{:d} duplications detected by Treerecs, extracting {:d} groups of at least {} orthologs.".format(len(dupl), 
																												  nDuplSign, 
																												  nbSp))

	return lOut
	
	
def runTreerecs(tree, sptree, o):
	"""
	Procedure which launch treerecs.

	@param1 tree: Path to the tree file
	@param2 sptree: Path to the species tree
	@param3 o: Output directory
	"""
	recTree = o+tree.split("/")[-1].replace(".txt","_recs.nhx")
	val = "treerecs -g {:s} -s {:s} -o {:s} -f -t 0.8 -O NHX:svg".format(tree, 
																	     sptree, 
																	     o)
	
	AnalysisFunc.cmd(val, False)
	
	return recTree


def treeTreatment(data, dAlnTree, nbSp, logger):
	"""
	Procedure which execute all functions for the tree step.

	@param1 data: basicData object
	@param2 logger: Logging object
	"""
	#try:
	dAlnTree2 = {}
	lFastaFile = []
	for aln, tree in dAlnTree.items():
		filtTree = filterTree(tree, 
							  data.sptree, 
							  data.cor)
		logger.info("Cleaning the tree")

		logger.info("Running Treerecs")
		recTree = runTreerecs(filtTree, 
							  data.sptree, 
							  data.o)
		#setattr(data, "recTree", recTree)

		lFastaFile += treeParsing(data.ORFs, 
								  recTree, 
								  nbSp, 
								  data.o, 
								  logger)
	
	setattr(data, "duplication", lFastaFile)
		
	if len(lFastaFile) > 0:
		for orthoGp in data.duplication:
			aln = AnalysisFunc.runPrank(orthoGp, 
									    data.geneName, 
									    data.o)
			tree = AnalysisFunc.runPhyML(aln, data.o)
			dAlnTree2[aln] = tree+"_phyml_tree.txt"

	dAlnTree.update(dAlnTree2)

	return dAlnTree

	#except Exception:
	#	logger.info("Treerecs uncountered an unexpected error, skipping.")
	#	return dicoAT


#######=================================================================================================================

