from blast_test import parseBlast
import logging
import pickle
import sys

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
	with open(sys.argv[1], 'rb') as fichier:
		data = pickle.load(fichier)
	
	accnFile = makeAccnsFile(data["lBlastRes"], sys.argv[2])
	data[accnFile] = accnFile

	with open(sys.argv[3],'wb') as fichier_data:
		pickle.dump(data,fichier_data,pickle.HIGHEST_PROTOCOL)