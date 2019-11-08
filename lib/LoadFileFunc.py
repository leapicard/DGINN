import BlastFunc, DataObject, ExtractFunc, TreeFunc, AnalysisFunc, FormatFunc
import os, sys, logging
from Bio import SeqIO

"""
This file pools the necessary functions to enter in the pipeline at a specific step.
"""

def checkPath(path, attrNames):
	"""
	Check if the path exist

	@param1 path: path's file 
	@param2 attrNames: Name of the attribute
	"""
	if not os.path.exists(path):
		print("The file with the {:s} information does not exist.".format(attrNames))
		sys.exit()

def baseNameInit(baseName, queryFile, aln, logger):
	"""
	Initialization of the attribut basename

	@param1 baseName: String
	@param2 queryFile: Path
	@param3 aln: Path
	@param4 logger: An object logging
	@return baseName: String
	"""
	
	if baseName == "":
		if queryFile != "":
			baseName = queryFile.split(".")[0].split("/")[-1]
		elif aln != "":
			baseName = aln.split(".")[0].split("/")[-1]
		else:
			logger.info("Basename can't be initialized.")
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
	corsg = TreeFunc.assocFile(sptree, filePath, o)

	path = TreeFunc.supData(filePath, corsg, o)

	return path, corsg

def parseLFile(path):
	"""
	Function which return a list of path from a file

	@param path: Path
	@return list of path
	"""
	with open(path, "r") as listFile:
		listFile = listFile.readlines()

	return [i.strip("\n") for i in listFile]

##==============================================================================================================================

def accnEntry(Data):
	"""
	Function handling start of the pipeline at the Extract step.

	@param Data: basicData object 
	@return data: basicData object
	"""
	logger = logging.getLogger("main")
	
	if FormatFunc.isBlastRes(Data.queryFile):
		Data.blastRes = Data.queryFile
		Data.lBlastRes = BlastFunc.parseBlast(Data.blastRes)
		Data.baseName = baseNameInit(Data.baseName, 
									 Data.queryFile, 
									 Data.aln, 
									 Data.logger)
	else:
		logger.info("The provided file is not a tabular output of Blast+, exiting DGINN.")
		sys.exit()
	
	return Data

def getSeqEntry(Data, treeOption):
	"""
	Function handling start of the pipeline at the Fasta step.

	@param Data: basicData object 
	@return data: basicData object
	"""
	logger = logging.getLogger("main")
	
	if FormatFunc.isAccns(Data.queryFile):
		Data.accnFile = Data.queryFile
		Data.baseName = baseNameInit(Data.baseName, 
									 Data.queryFile, 
									 Data.accnFile, 
									 Data.logger)

		Data.lBlastRes = [ i.strip("\n") for i in open(Data.accnFile, "r").readlines() ]
	else:
		logger.info("Provided file is not a list of NCBI accessions, terminating DGINN.")
		sys.exit()
		
	return Data


def orfEntry(Data, treeOption):
	"""
	Function handling start of the pipeline at the orf step.

	@param1 Data: basicData object
	@param2 treeOption: Boolean
	@return data: basicData object
	"""
	logger = logging.getLogger("main")
	
	if FormatFunc.isFasta(Data.queryFile):
		Data.seqFile = Data.queryFile
		Data.baseName = baseNameInit(Data.baseName, 
									 Data.queryFile, 
									 Data.aln, 
									 Data.logger)
	else:
		logger.info("The provided file is not a fasta of nucleotide sequences, exiting DGINN.")
		sys.exit()
	
	return Data

def prankEntry(Data, treeOption):
	"""
	Function handling start of the pipeline at the Prank step.

	@param1 Data: basicData object
	@param2 treeOption: Boolean
	@return data: basicData object
	"""
	if FormatFunc.isFasta(Data.queryFile):
		Data.ORFs = Data.queryFile
		Data.baseName = baseNameInit(Data.baseName, 
									 Data.queryFile, 
									 Data.ORFs, 
									 Data.logger)
		
		with open(Data.ORFs) as orf:
			Data.geneName = orf.readline().split("_")[1]
		
	else:
		logger.info("Provided file is not a fasta of sequences, terminating DGINN.")
		sys.exit()
		
	return Data

def phymlRecEntry(Data, logger):
	"""
	Function handling start of the pipeline at the phyml step.

	@param1 Data: basicData object
	@param2 treeOption: Boolean
	@param3 logger: An object logging
	@return Data: basicData object
	"""
	if FormatFunc.isAln(Data.queryFile):
		Data.aln = Data.queryFile
		Data.ORFs = Data.queryFile
		Data.baseName = baseNameInit(Data.baseName, 
									 Data.queryFile, 
									 Data.aln, 
									 Data.logger)
		
		with open(Data.aln) as orf:
			Data.geneName = orf.readline().split("_")[1]
		
	else:
		logger.info("Provided file is not a multiple sequence alignment, terminating DGINN.")
		
	return Data

def spTreeCheck(Data, firstStep, treeOption):
	if treeOption:
		checkPath(Data.sptree, "species's tree")
		
		if not hasattr(Data, 'cor'):
			if firstStep == "orf":
				Data.seqFile, corSG = filterData(Data.sptree, 
												 Data.seqFile, 
												 Data.o)
			elif firstStep == "prank":
				Data.ORFs, corSG = filterData(Data.sptree, 
											  Data.ORFs, 
											  Data.o)
			elif firstStep == "phyml":
				Data.aln, corSG = filterData(Data.sptree, 
											 Data.aln, 
											 Data.o)
			elif firstStep == "duplication":
				Data.aln, corSG = filterData(Data.sptree, 
											 Data.aln, 
											 Data.o)
			setattr(Data, "cor", corSG)

### rework this
def duplPSEntry(Data, logger):
	"""
	Function handling start of the pipeline at the tree step.

	@param Data: basicData object
	@param2 logger: Logging object
	@return Data: basicData object
	"""
	dico = {}
	if Data.aln != "" and Data.tree != "":
		logger.info("Alignement file: "+Data.aln)
		logger.info("Gene Tree file: "+Data.tree)
		dico[Data.aln] = Data.tree
		Data.ORFs = Data.aln
		Data.baseName = baseNameInit(Data.baseName, Data.queryFile, Data.aln, Data.logger)
	else:
		logger.info("Alignment and/or gene tree file have not been provided.")
		sys.exit()

	return Data, dico

def gardEntry(Data, parameters, logger):
	"""
	Function handling start of the pipeline at the GARD step.

	@param1 Data: basicData object
	@param2 parameters: Dico of parameters
	@param3 logger: Logging object
	@return Data: basicData object
	"""
	dico = {}
	Data = communFuncEntry(Data, [], 1)
	#AccessionFunc.createGeneDir(Data.o, Data.aln.split('/')[-1].split(".")[0])
	if Data.tree != "" and Data.aln != "":
		dico[Data.aln] = Data.tree
	else:
		logger.info("You didn't precise the alnfile or the treefile.")
		sys.exit()
	return Data, dico


def pspEntry(Data, parameters, logger):
	Data, dico = LoadFileFunc.duplPSEntry(Data, logger)
	
	Data.alnFormat = parameters["alnformat"].title()

	return Data, dico

##==============================================================================================================================
