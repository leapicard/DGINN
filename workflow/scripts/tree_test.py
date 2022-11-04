import sys
import logging, os
import ete3
import fastares_test
from Bio import SeqIO
from collections import defaultdict

"""
File which countain all functions about treerecs and tree treatement.
"""
def getLeaves(path):
	"""
	Open a newick file and return a list of species

	@param path: path of a tree file
	@return lTreeData: list of species
	"""
	tree = ete3.Tree(path)
	lGene = tree.get_leaf_names()

	return lGene

#=========================================================================================================================

def assocFile(sptree, path, dirName):
	"""
	Create a file which contain the species of each genes

	@param1 sptree: Species tree
	@param2 path: Path of a fasta file
	@param3 dirName: Name for a new directory
	"""        

	lGeneId = []
	ff=open(path, "r")
	for accn in SeqIO.parse(ff, "fasta"):
		lGeneId.append(accn.id)
	ff.close()
	lGeneId.sort()
	        
	lTreeData = getLeaves(sptree)
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
		for k,v in dSp2Gen.items():
			cor.write(k+"\t"+v+"\n")
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
		lGene = [ i.split("\t")[0] for i in lcorSG ]
		corSG.close()
	newDico = {}
	for accn in SeqIO.parse(open(filePath, "r"), "fasta"):
		if accn.id in lGene:
			newDico[accn.id] = accn.seq

	out = dirName+filePath.replace(".fasta", "_filtered.fasta").split("/")[-1]
	with open(out, "w") as newVer:
		newVer.write(fastares_test.dict2fasta(newDico))
		newVer.close()
	return out

def buildSpeciesTree(taxon, gfaln):
	"""
	Build a species tree from the ncbi taxonomy and species in the gene tree,
	and write it in a file.

	@param1 taxon: taxon under which ncbi species are considered.
	@param2 gftree: path of the alignment
	@return the species tree file path.
	"""

	ncbi = ete3.NCBITaxa(dbfile="/opt/DGINN/taxa.sqlite")

	sptree = ncbi.get_descendant_taxa(taxon, collapse_subspecies=True, return_tree=True)
	spleaves = [ncbi.get_taxid_translator([x]) for x in sptree.get_leaf_names()]

	accns = list(SeqIO.parse(gfaln,'fasta'))
	gLeavesSp = [name.id.split("_")[0] for name in accns]

	## link Taxids & species names abbreviations
	leavesNames = {}
	for x in spleaves:
		for k, v in x.items():
			tax = v.replace(".", "").replace("-", "").split(" ")
			newTax = tax[0][:3].lower()+"".join([ i[:3].title() for i in tax[1:]])
			leavesNames[k] = newTax

	# List of tax ids of species in gene tree
	lTaxids = []
	for x in gLeavesSp:
		lTaxids += [str(k) for k, v in leavesNames.items() if v == x[0:6]]

	lTaxids = list(set(lTaxids))

	# restriction of species tree to these taxons
	sptree.prune(lTaxids)

	# back to correct leaves names
	for x in sptree:
		x.name = leavesNames[int(x.name)]

	spTreeFile = "/".join(gfaln.split("/")[:-1]+["species_tree.tree"])
	sptree.write(format=9, outfile=spTreeFile)

	return spTreeFile

  
def treeCheck(treePath, alnf, optionTree):
	"""
	Check if the tree isn't corruped

	@param1 treePath: tree's path
	@param2 alnf: alignment's path
	@param3 optionTree: Boolean (treerecs option)
	@return1 treePath: tree's path
	@return2 optionTree: Boolean (treerecs option)
	"""
	logger = logging.getLogger("main")
	if treePath != "":
		if not os.path.exists(treePath):
			try:
				treePath = buildSpeciesTree(treePath, alnf)
			except:
				logger.warning("The path for the species tree doesn't exist.")
				logger.warning("Duplication option will be set to False.")
				optionTree = False
				treePath = ""
		else:
			if not ete3.Tree(treePath):
				logger.warning("The species tree is corrupted.")
				logger.warning("Duplication option will be set to False.")
				optionTree = False
				treePath = ""

	return treePath, optionTree