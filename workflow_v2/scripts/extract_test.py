from blast_test import parseBlast
import logging
import json
import sys
import loadfile_test

"""This file pools the necessary functions to treat the input file of genes and their CCDS accessions."""

def makeAccnsFile(lBlastRes, output_file):
	"""
	Function creating a file with only CCDS accessions.

	@param1 lBlastRes: List of accessions
	@param2 queryName: ID of the blast query/reference sequence
	@param3 outDir: Output directory
	@return outBlastn: Path to the accessions file
	"""

	logger = logging.getLogger("main.accessions")
	geneAllAccns = "\n".join(set(lBlastRes))
	
	logger.debug("All Accns: " + ",".join(set(lBlastRes)))

	with open(output_file, "w") as out:
		out.write(geneAllAccns)
		out.close()	
	logger.info("Wrote all accessions to {:s}".format(output_file))
	return(output_file)

if __name__ == "__main__" :

	with open(sys.argv[1], 'r') as config_in:
		config_dict = json.load(config_in)
	
	if config_dict["parameters"]["step"] == "accessions" :
		config_dict["data"] = loadfile_test.accnEntry(config_dict["data"])

	accnFile = makeAccnsFile(config_dict["data"]["lBlastRes"], sys.argv[3])		#!#
	config_dict["data"]["accnFile"] = accnFile

	with open(sys.argv[1],'w') as config_out:
		json.dump(config_dict, config_out)
		