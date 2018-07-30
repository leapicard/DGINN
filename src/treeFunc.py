import logging, os
import ete3
import FastaResFunc, AnalysisFunc, AccessionFunc
from Bio import SeqIO

"""
File which countain all functions about treerecs and tree treatement.
"""

#######Treerecs=========================================================================================================
#=========================================================================================================================
def speciesTree(path):
	"""
	Open a newick file and return a list of species

	@param path: path of a tree file
	@return lTreeData: list of species
	"""
	with open(path, "r") as tree:
		treeData = tree.readlines()[0]
		lGene = getSp(treeData)

	return lGene

def getSp(string):
	"""
	Function which take a tree to extract and extract all word (suppression of (;): and ,) and return a list

	@param string: A tree (newick format)
	@return ID: List of words 
	"""
	lTreeData = string.replace("(", ",").replace(")", ",").replace(":", ",").replace(";", ",").split(",")
	ID = [ i for i in lTreeData[:-1] if i != "" and i[1] != "." ]

	return ID

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
			logger.info("The path for the species tree doesn't exist.")
			logger.info("Treerecs option will be set to False.")
			optionTree = False
		else:
			with open(treePath, "r") as sptree:
				sptree = sptree.readline()
				nbOpening = sptree.count("(")
				nbEnding = sptree.count(")")
				if nbOpening != nbEnding:
					logger.info("The species's tree is corruped...")
					logger.info("Treerecs option will be set to False.")
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
		
	lTreeData = speciesTree(sptree)

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

	out = dirName+filePath.replace(".fasta", "_filtred.fasta").split("/")[-1]
	with open(out, "w") as newVer:
		newVer.write(FastaResFunc.dict2fasta(newDico))

	return out

#=========================================================================================================================

def treeFiltred(tree, spTree, cor):
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
		dCor[i[0]] = i[1].strip("\n")
		dCorInv[i[1].strip("\n")] = i[0]

	g1Genes = []
	for z in g1:
		if dCor[z] in g2:
			g1Genes.append(z)

	if len(g1Genes) < len(g1)/2:
		print(len(g1Genes))
		print(len(g1)/2)
		exit()

	print(g1Genes)
	t1.prune(set(g1Genes))
	out = tree.replace(".tree", "_filtred.tree")
	t1.write(out)
	return out

#=========================================================================================================================
def supDuplic(cor, tree):
	"""
	Deleting of genes whose species is the same and are "neighbour"

	@param1 cor: Path to the file with correspondence beetween gene and species
	@param2 tree: Path to a tree file
	@return out: Path
	"""
	dCor = {}
	for line in open(cor, "r"):
		temp = line.split("\t")
		dCor[temp[0]] = temp[1].strip("/")

	with open(tree, "r") as t:
		t = t.readlines()[0].strip("\n")

	tMod = t
	temp = t.replace("(",")").split(")")
	lSup = []
	#We split by "(" and ")" and we search which tuple have a lengh equal to 2. If there're 2 gene with the same species, we delete one of the two in the tree
	for j in temp:
		lValor = j.split(",")
		if len(lValor) == 2:
			if "" not in lValor:
				gene1 = lValor[0].split(":")[0]
				gene2 = lValor[1].split(":")[0]
				if gene1 in dCor and gene2 in dCor:
					if dCor[gene1] == dCor[gene2]:
						lSup.append(lValor[1])

	for i in lSup:
		tMod = tMod.replace(","+i, "")

	out = tree.replace(".txt", "_modified.txt")

	with open(out, "w") as treeMod:
		treeMod.write(tMod)

	return out

def treeParsing(ORFs, tree, geneName, o, logger):
	"""
	Function which parse gene data in many group according to duplication in the reconciliated tree

	@param1 ORFs: Path to the ORFs file
	@param2 tree: Path to a tree file
	@param3 geneName: Gene name
	@param4 o: Output directory
	@param5 logger: Logging object
	@return lOut: List of path (fasta files)
	"""
	with open(o+tree.split("/")[-1].replace(".txt","_reconciliation.newick"), "r") as tree:
		reconTree = tree.readlines()[1]

	#Looking for the word duplication and save the index
	lIndex = []
	index = 0
	temp = -1
	while index != -1:
		index = reconTree.find("duplication", temp + 1)
		temp = index
		if index != -1:
			lIndex.append(index)

	logger.info("Treerecs has identified {:s} duplications in the alignement for {:s}.".format(str(len(lIndex)), geneName))

	#Save subtree of each duplication
	lGroup = []
	for i in lIndex:
		j = i-2
		begin = 0
		end = 1
		while j != 0:
			if reconTree[j] == ")":
				end += 1
			elif reconTree[j] == "(":
				begin += 1
			j -= 1
			if begin-end == 0:
				break
		lGroup.append(reconTree[j+2:i-1])

	#Divsion of each subTree in two group of genes
	number = 1
	lOut = []
	counter = 0
	for j in lGroup:
		#While the difference beetween the number of"(" and ")" is difference to 0, we continue to read the tree
		count = 0
		index2 = 1
		if "(" in j and ")" in j:
			while count != -42:
				if j[index2] == "(":
					count += 1
				elif j[index2] == ")":
					count -= 1
				if count == 0:
					count = -42
				index2 += 1

			list1 = getSp(j[:index2])
			list2 = getSp(j[index2+1:])
			list1 = [i for i in list1 if "_loss" not in i and "duplication" not in i]
			list2 = [i for i in list2 if "_loss" not in i and "duplication" not in i]
		else:
			list1 = []
			list2 = []

		dFasta1 = {}
		dFasta2 = {}
		if len(list1) >8 and len(list2) > 8:
			counter += 1
			for gene in SeqIO.parse(open(ORFs),'fasta'):
				if gene.id in list1:
					dFasta1[gene.id] = gene.seq
				elif gene.id in list2:
					dFasta2[gene.id] = gene.seq

			#Creation of direcotory for each group and duplication, and safeguard
			o1 = AccessionFunc.createGeneDir(o, "D"+str(number)+"_G1")
			o2 = AccessionFunc.createGeneDir(o, "D"+str(number)+"_G2")
			out = o1+ORFs.split("/")[-1].split(".")[0]+"_D"+str(number)+"_G1.fasta"
			out2 = o2+ORFs.split("/")[-1].split(".")[0]+"_D"+str(number)+"_G2.fasta"
			with open(out, "w") as fastaGr1:
				fastaGr1.write(FastaResFunc.dict2fasta(dFasta1))

			with open(out2, "w") as fastaGr2:
				fastaGr2.write(FastaResFunc.dict2fasta(dFasta2))
			number += 1
			lOut.append(out)
			lOut.append(out2)

		elif len(list1) >8 or len(list2) >8:
			if len(list1) >8:
				lg = list1
			else:
				lg = list2

			for gene in SeqIO.parse(open(ORFs),'fasta'):
				if gene.id in lg:
					dFasta1[gene.id] = gene.seq
			o1 = AccessionFunc.createGeneDir(o, "D"+str(number))
			out = o1+ORFs.split("/")[-1].split(".")[0]+"_D"+str(number)+".fasta"
			with open(out, "w") as fasta:
				fasta.write(FastaResFunc.dict2fasta(dFasta1))
			number += 1
			lOut.append(out)

	logger.info("{:d} duplications have been detected.".format(counter))

	return lOut

def alnsorting(reconTree, alnFile):
	"""
	Updating aln file by deleting genes whose aren't in the reconciliation tree

	@param1 reconTree: Path to the reconciliation tree
	@param2 alnFile: Path to the alignment file
	@return out: Path
	"""
	with open(reconTree, "r") as sptreerec:
		tree = sptreerec.readlines()[1]

	lSp = getSp(tree)
	newLSpGene = [i for i in lSp if "_loss" not in i and i != "duplication"]

	dicoAln = {}
	for gene in SeqIO.parse(open(alnFile),'fasta'):
		dicoAln[gene.id] = gene.seq

	newDicoAln = {}
	for gene in newLSpGene:
		newDicoAln[gene] = dicoAln[gene]

	out = alnFile.split(".")[0]+"_aligned_reconTree.fasta"

	with open(out, "w") as newFasta:
		newFasta.write(FastaResFunc.dict2fasta(newDicoAln))

	return out

def runTreerecs(tree, sptree, o):
	"""
	Procedure which lunch treerecs.

	@param1 tree: Path to the tree file
	@param2 sptree: Path to the species tree
	@param3 o: Output directory
	"""
	AnalysisFunc.cmd("Treerecs -g {:s} -s {:s} -R -o {:s}".format(tree, sptree, o+tree.split("/")[-1].replace(".txt","_reconciliation.newick")), False)

def runPrank2(duplication, logger):
	"""
	Function which lunch prank on each subgroup.

	@param1 duplication: List of path (fasta files)
	@param2 logger: Logging object
	@return lRes: List of path (alignments files)
	"""
	lRes = []
	for lFasta in duplication:
		outPrank = lFasta.replace(".fasta", "_prank")
		geneName = lFasta.split("/")[-1].split(".")[0]
		AnalysisFunc.runPrank(lFasta, geneName, outPrank)
		lRes.append(outPrank+".best.fas")

	return lRes

def treeTreatment(data, dicoAT, logger):
	"""
	Procedure which execute all functions for the tree step.

	@param1 data: basicData object
	@param2 logger: Logging object
	"""
	try:
		data.tree = treeFiltred(data.tree, data.sptree)
		logger.info("Cleaning the tree...")
		data.tree = supDuplic(data.cor, data.tree)

		logger.info("Beginning of treerecs works")
		runTreerecs(data.tree, data.sptree, data.o)

		lFastaFile = treeParsing(data.ORFs, data.tree, data.geneName, data.o, logger)
		setattr(data, "duplication", lFastaFile)

		data.aln = alnsorting(data.o+data.tree.split("/")[-1].replace(".txt","_reconciliation.newick"), data.aln)

		lResPrank = runPrank2(data.duplication, logger)
		
		for i in data.lResPrank:
			directory = "/".join(i.split("/")[:-1])
			tree = runPhyML(i, directory)
			data.lTreeAln.append(tree)
			dicoAT[i] = tree+"_phyml_tree.txt"

		return dicoAT

	except Exception:
		logger.info("Treerecs uncountered an unexpected error, skipping.")
		return dicoAT


#######=================================================================================================================

