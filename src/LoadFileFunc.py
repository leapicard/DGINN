import Blast, G2P_Object, AccessionFunc, treeFunc, AnalysisFunc
import os
from multiprocessing import Pool
from Bio import SeqIO

"""
This file pools the necessary functions to enter in the pipeline at a specific step.
"""

def securityEntry(lFiles, nbFiles):
	"""
	Verify if the number of files.

	@param1 lFiles: list of paths
	@param2 nbFiles: number of files needed

	"""
	if len(lFiles) != nbFiles:
		print("There're not enough files or too much to enter at this step !")
		exit()

def checkPath(path, attrNames):
	"""
	Check if the path exist

	@param1 path: path's file 
	@param2 attrNames: Name of the attribut
	"""
	if not os.path.exists(path):
		print("The File with the {:s} informations does not exist.".format(attrNames))
		exit()

def baseNameInit(baseName, CCDSFile, aln, logger):
	"""
	Initialization of the attribut basename

	@param1 baseName: String
	@param2 CCDSFile: Path
	@param3 aln: Path
	@param4 logger: An object logging
	@return baseName: String
	"""
	if baseName == "":
		if CCDSFile != "":
			baseName = CCDSFile.split(".")[0].split("/")[-1]
		elif aln != "":
			baseName = aln.split(".")[0].split("/")[-1]
		else:
			logger.info("We can't initialize the basename")
			exit()
	return baseName

def filtreData(sptree, filePath, o):
	"""
	Function which execute functions to delete genes which aren't in the species tree.

	@param1 sptree: Path of a newick file
	@param2 filePath: Path of a Fasta file
	@param3 o: Path of a directory
	@return path: Path of a file
	"""
	corsg = treeFunc.assocFile(sptree, filePath, o)

	path = treeFunc.supData(filePath, corsg, o)

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
def communFuncEntry(Data, attrNames, nb):
	"""
	Commun code for function in LoadFileFunc.py.

	@param1 Data: basicData object
	@param2 attrNames: Name of the new attribut
	@param3 nb: Number of files
	@return data: basicData object
	"""
	files = Data.CCDSFile.split(",")
	securityEntry(files, nb)
	Data.CCDSFile = files[0].strip("\n")
	checkPath(Data.CCDSFile, "fasta")
	Data.setGenAttr()

	setattr(Data, "geneDir", Data.o+Data.geneName.split("|")[1]+"/")
	os.makedirs(Data.geneDir)

	for i in range(len(attrNames)):
		checkPath(files[i+1].strip("\n"), attrNames[i])
		setattr(Data, attrNames[i], files[i+1].strip("\n"))

	Data.baseName = baseNameInit(Data.baseName, Data.CCDSFile, Data.aln, Data.logger)

	return Data

def accessionEntry(Data):
	"""
	Function handling start of the pipeline at the Accession step.

	@param Data: basicData object 
	@return data: basicData object
	"""

	files = Data.CCDSFile.split(",")
	securityEntry(files, 2)
	dId = {}
	Data.CCDSFile = files[0].strip("\n")
	checkPath(Data.CCDSFile, "fasta")
	Data.setGenAttr()
	Data.blastRes = files[1].strip("\n")
	checkPath(Data.blastRes, "Blast' resultats")
	Data.lBlastRes = Blast.parseBlast(Data)
	Data.baseName =  baseNameInit(Data.basename, Data.CCDSFile, Data.aln, Data.logger)

	return Data

def fastaEntry(Data):
	"""
	Function handling start of the pipeline at the Fasta step.

	@param Data: basicData object 
	@return data: basicData object
	"""
	Data = communFuncEntry(Data, ["accnFile"], 2)

	Data.lBlastRes = [ i.strip("\n") for i in open(Data.accnFile, "r").readlines() ]

	return Data


def orfEntry(Data, treeOption):
	"""
	Function handling start of the pipeline at the orf step.

	@param1 Data: basicData object
	@param2 treeOption: Boolean
	@return data: basicData object
	"""
	Data = communFuncEntry(Data, ["catFile"], 2)
	if treeOption == "True":
		checkPath(Data.sptree, "species's tree")
		Data.catFile, corSG = filtreData(Data.sptree, Data.catFile, Data.o)
		setattr(Data, "cor", corSG)
	return Data

def prankEntry(Data, treeOption):
	"""
	Function handling start of the pipeline at the Prank (currently Mafft) step.

	@param1 Data: basicData object
	@param2 treeOption: Boolean
	@return data: basicData object
	"""
	Data = communFuncEntry(Data, ["catFile", "ORFs"], 3)
	if treeOption == "True":
		checkPath(Data.sptree, "species's tree")
		Data.ORFs, corSG = filtreData(Data.sptree, Data.ORFs, Data.o)
		setattr(Data, "cor", corSG)
	return Data

def phymlEntry(Data, treeOption, logger):
	"""
	Function handling start of the pipeline at the phyml step.

	@param1 Data: basicData object
	@param2 treeOption: Boolean
	@param3 logger: An object logging
	@return Data: basicData object
	"""
	if Data.aln != "":
		logger.info("Alignement file: "+Data.aln)
		Data = communFuncEntry(Data, ["ORFs"], 2)
	else:
		logger.info("You need to precise the alnfile")
		exit()
	if treeOption == "True":
		checkPath(Data.sptree, "species's tree")		
		Data.aln, corSG = filtreData(Data.sptree, Data.aln, Data.o)
		setattr(Data, "cor", corSG)
		Data.tree = AnalysisFunc.runPhyML(Data.aln, "/".join(Data.aln.split("/")[:-1]))+"_phyml_tree.txt"
	return Data


def treeEntry(Data, logger):
	"""
	Function handling start of the pipeline at the tree step.

	@param Data: basicData object
	@param2 logger: Logging object
	@return Data: basicData object
	"""
	dico = {}
	if Data.aln != "" and Data.tree != "":
		logger.info("Alignement file: "+Data.aln)
		logger.info("Gene's Tree file: "+Data.tree)
		Data = communFuncEntry(Data, ["ORFs"], 2)
	else:
		logger.info("You didn't precise the alnfile or the treefile.")
		exit()

	checkPath(Data.sptree, "species's tree")		
	Data.aln, corSG = filtreData(Data.sptree, Data.aln, Data.o)
	setattr(Data, "cor", corSG)
	dico[Data.aln] = Data.tree
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
	AccessionFunc.createGeneDir(Data.o, Data.aln.split('/')[-1].split(".")[0])
	if Data.tree != "" and Data.aln != "":
		dico[Data.aln] = Data.tree
	else:
		logger.info("You didn't precise the alnfile or the treefile.")
		exit()
	return Data, dico


def pspEntry(Data, parameters, logger):
	dico = {}
	logger.info("Alignement file: "+Data.aln)
	logger.info("Gene's Tree file: "+Data.tree)
	dico[Data.aln] = Data.tree
	Data.baseName = baseNameInit(Data.baseName, Data.CCDSFile, Data.aln, Data.logger)
	AccessionFunc.createGeneDir(Data.o, Data.aln.split('/')[-1].split(".")[0])
	Data.alnFormat = parameters["alnformat"].title()

	return Data, dico

##==============================================================================================================================
