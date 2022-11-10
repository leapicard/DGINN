import tree_test, FormatFunc, blast_test
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

def formatcheck_isFasta(targetFile):
	"""
	Check if provided file is a simple fasta and contain multiple sequences.
	
	@param targetFile: file to be tested
	@return boolFile: boolean of answer
	"""
	
	accns = list(SeqIO.parse(targetFile,'fasta'))
	boolFile = any(accns)
	
	if len(accns) < 2:
		boolFile = False
	
	return(boolFile)

def accnEntry(data_dict):
	"""
	Function handling start of the pipeline at the Extract step.

	@param Data: basicData object 
	@return data: basicData object
	"""
	if FormatFunc.isBlastRes(data_dict["queryFile"]):
		data_dict["blastRes"] = data_dict["queryFile"]
		data_dict["lBlastRes"] = blast_test.parseBlast(data_dict["blastRes"])
		data_dict["baseName"] = baseNameInit(data_dict["baseName"], 
							data_dict["queryFile"], 
							data_dict["aln"],
							"accessions")
	else:
		logger = logging.getLogger("main.accessions")
		logger.error("The provided file is not a tabular output of Blast+, exiting DGINN.")
		sys.exit()

	return data_dict

# def orfEntry(data_dict):
# 	"""
# 	Function handling start of the pipeline at the orf step.

# 	@param1 Data: basicData object
# 	@return data: basicData object
# 	"""
# # Condition checking if queryFile is a fasta file r
# 	if formatcheck_isFasta(data_dict["queryFile"]):
# 		data_dict["seqFile"] = data_dict["queryFile"]
# 		data_dict["baseName"] = baseNameInit(data_dict["baseName"], 
# 					     data_dict["queryFile"], 
# 					     data_dict["aln"],
#                         "orf")
# 	else:
# 		logger = logging.getLogger("main.orf")
# 		logger.error("The provided file is not a fasta of nucleotide sequences, exiting DGINN.")
# 		sys.exit()
	
# 	return data_dict