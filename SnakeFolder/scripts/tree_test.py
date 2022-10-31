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
