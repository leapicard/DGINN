import tree_test
import os, sys, logging
from Bio import SeqIO

def baseNameInit(baseName, queryFile, aln, step = ""):
	"""
	Initialization of the attribut basename

	@param1 baseName: String
	@param2 queryFile: Path
	@param3 aln: Path
        @param4
	@return baseName: String
	"""
	
	logger = logging.getLogger(".".join(["main",step]))
	if baseName == "":
		if queryFile != "":
			baseName = queryFile.split(".")[0].split("/")[-1]
		elif aln != "":
			baseName = aln.split(".")[0].split("/")[-1]
		else:
			logger.error("Basename can't be initialized.")
			sys.exit()
	return baseName
	
def filterData(sptree, filePath, o):
	"""
	Function which execute functions to delete genes which aren't in the species tree.

	@param1 sptree: Path of a newick file
	@param2 filePath: Path of a Fasta file
	@param3 o: Path of a directory
	@return path: Path of a file
	"""
	
	corsg = tree_test.assocFile(sptree, filePath, o)

	path = tree_test.supData(filePath, corsg, o)

	return path, corsg

