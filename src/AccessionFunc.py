import logging, os
import G2P_Object
from multiprocessing import Pool
from collections import defaultdict, OrderedDict

"""This file pools the necessary functions to treat the input file of genes and their CCDS accessions."""

def createGeneDir(o, Name):
	"""
	Function creating the output directory where results will be filed.

	@param1 o: Output directory
	@param2 Name: Directory name
	@return geneDir: Path to the output directory
	"""

	# create subdirectory for gene results
	logger = logging.getLogger("main")
	geneDir = o+Name+"/"
	if not os.path.exists(geneDir):
		os.makedirs(geneDir)
	
	logger.info("Results for {:s} will be output to {:s}".format(Name, geneDir))

	return(geneDir)

def makeAccnsFile(lBlastRes, geneName, geneDir):
	"""
	Function creating a file with only CCDS accessions.

	@param1 lBlastRes: List of accessions
	@param2 geneName: Gene name
	@param3 geneDir: gene Directory
	@return outBlastn: Path to the accessions file
	"""

	# write all accessions to new file		
	outBlastn = geneDir+geneName.split("|")[1]+"_accns.txt"
	geneAllAccns = ""
	logger = logging.getLogger("main")

	geneAllAccns = "\n".join(lBlastRes)
	
	logger.debug(geneAllAccns)
		
	with open(outBlastn, "w") as out:
		out.write(geneAllAccns)

	out.close()
	
	logger.info("{:s}: wrote all accessions to {:s}".format(geneName.split("|")[1], outBlastn))
	
	return(outBlastn)


def treatAccns(Data, logger):
	"""
	Process handling every accession-related step of the pipeline.

	@param1 Data: basicdata object
	@param2 logger: Logging object
	"""

	geneDir = createGeneDir(Data.o, Data.geneName.split("|")[1])
	setattr(Data, "geneDir", geneDir)
	
	logger.info("Wrote accessions to file")
	accnFile = makeAccnsFile(Data.lBlastRes, Data.geneName, Data.geneDir)
	setattr(Data,"accnFile", accnFile)
	
