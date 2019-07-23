import logging, os
import G2P_Object
from collections import defaultdict, OrderedDict

"""This file pools the necessary functions to treat the input file of genes and their CCDS accessions."""


def makeAccnsFile(lBlastRes, geneName, outDir):
	"""
	Function creating a file with only CCDS accessions.

	@param1 lBlastRes: List of accessions
	@param2 geneName: Gene name
	@param3 geneDir: gene Directory
	@return outBlastn: Path to the accessions file
	"""

	# write all accessions to new file		
	outBlastn = outDir+geneName.split("|")[1]+"_accns.txt"
	geneAllAccns = ""
	logger = logging.getLogger("main")

	geneAllAccns = "\n".join(set(lBlastRes))
	
	logger.debug(geneAllAccns)
		
	with open(outBlastn, "w") as out:
		out.write(geneAllAccns)

	out.close()
	
	logger.info("Wrote all accessions to {:s}".format(outBlastn))
	
	return(outBlastn)


def treatAccns(Data, logger):
	"""
	Process handling all the accession-related steps of the pipeline.

	@param1 Data: basicdata object
	@param2 logger: Logging object
	"""
	
	accnFile = makeAccnsFile(Data.lBlastRes, Data.geneName, Data.o)
	setattr(Data,"accnFile", accnFile)
	
