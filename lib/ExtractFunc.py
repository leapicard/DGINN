#coding: utf8
import logging

"""This file pools the necessary functions to treat the input file of genes and their CCDS accessions."""


def makeAccnsFile(lBlastRes, queryName, outDir):
	"""
	Function creating a file with only CCDS accessions.

	@param1 lBlastRes: List of accessions
	@param2 queryName: ID of the blast query/reference sequence
	@param3 outDir: Output directory
	@return outBlastn: Path to the accessions file
	"""

	# write all accessions to new file
	outBlastn = outDir+queryName.split("_")[0:2][-1]+"_accns.txt"
	logger = logging.getLogger("main.accessions")
	geneAllAccns = "\n".join(set(lBlastRes))
	
	logger.debug(geneAllAccns)
		
	with open(outBlastn, "w") as out:
		out.write(geneAllAccns)
	
	logger.info("Wrote all accessions to {:s}".format(outBlastn))
	
	return(outBlastn)


def treatAccns(Data):
	"""
	Process handling all the accession-related steps of the pipeline.

	@param1 Data: basicdata object
	@param2 logger: Logging object
	"""
	
	accnFile = makeAccnsFile(Data.lBlastRes, Data.queryName, Data.o)
	setattr(Data,"accnFile", accnFile)
	
