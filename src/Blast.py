import G2P_Object, AnalysisFunc
import logging, os, shlex
from collections import defaultdict, OrderedDict
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SearchIO

"""This file pools the necessary functions to run Blast and parse its results."""

def blast(CCDSFile, db, evalue, percId, cov, apiKey, remote, query, logger):
	"""
	Function running Blast on the provided gene list.
	
	@param1 CCDSFile: Path
	@param2 db: Database name
	@param3 evalue: Blast e-value
	@param4 percId: Blast percentage of identity
	@param5 cov: Blast coverage
	@param6 apiKey: Key for the API of NCBI
	@param7 remote: Boolean (online database or not)
	@param8 query: Strings which will select the database online
	@param9 logger: Logging object
	@return blastRes: Path to Blast results file
	"""

	# Blast results file
	blastRes = CCDSFile.split(".")[0]+"_blastres.tsv"

	logger.info("Blasting sequence...")

	if remote == "True":
		
		sequences = open(CCDSFile).read()	
		resultHandle = NCBIWWW.qblast("blastn", db, sequences, entrez_query= query, url_base='https://blast.ncbi.nlm.nih.gov/Blast.cgi?api_key={:s}'.format(apiKey),hitlist_size=10000, perc_ident= percId, threshold = evalue)
		
		f = open(blastRes, "w")
		for record in NCBIXML.parse(resultHandle):
			f.write("# "+record.application+" "+record.version+"\n# Query: "+record.query+"\n# Database: "+record.database+"\n# Fields: subject id, ?\n# "+str(len(record.alignments))+" hits found\n")
			for alignment in record.alignments:
				for hsp in alignment.hsps:
					f.write(alignment.title+"\n")
		f.close()

	else:
		AnalysisFunc.cmd("blastn -task blastn -db {:s} -query {:s} -out {:s} -evalue {:f} -perc_identity {:f} -qcov_hsp_perc {:f} -outfmt \"7 sseqid\" -max_target_seqs 10000 -num_threads {:d}".format(db, CCDSFile, blastRes, evalue, percId, cov, 1), False)
		logger.debug("Blast command: {:s}".format(cmd))


	if os.path.exists(blastRes) or os.stat(blastRes).st_size > 0:								
		logger.info("Blast results written to: {:s}".format(blastRes))
	else:
		logger.error("Blast didn't run, stopping script.")
		exit()

	return(blastRes)

def parseBlast(blastRes, logger):
	"""
	Function parsing accessions of Blast results per gene.

	@param1 blastRes: Path
	@param2 logger: Logging object
	@return lGeneRes: list of all the accessions of homologous genes found through Blast
	"""

	# split blast results by gene using first line of comments as a separator
	with open(blastRes, "r") as blast:
		separator = blast.readline().strip()
		listBlastRes = blast.read()
		logger.debug(listBlastRes)
		
		lGeneRes = defaultdict(list)
		listGeneRes = list(filter(bool, listBlastRes.split("\n")))
		geneName = listGeneRes[0].split("|")[1]
		if listGeneRes[-1][0] == "#":
			listGeneRes = listGeneRes[:-1]

		lGeneRes = [ i.split("\t")[0].split("|")[-2]for i in listGeneRes[5:] ]

	logger.debug(lGeneRes)
	logger.info("Separated blast results per gene")

	return(lGeneRes)

def treatBlast(data, evalue, percId, cov, apiKey, remote, query):
	"""
	Process creating a list of gene objects (one for each gene of the input list).

	@param1 data: A basicData object
	@param2 evalue: Blast e-value
	@param3 percId: Blast percentage of identity
	@param4 cov: Blast coverage
	@param5 apiKey: Key for the API of NCBI
	@param6 remote: Boolean (online database or not)
	@param7 query: Strings which will select the database online
	@return data: Liste of gene objects
	"""

	data.blastRes = blast(data.CCDSFile, data.db, evalue, percId, cov, apiKey, remote, query, data.logger)
	data.lBlastRes = parseBlast(data.blastRes, data.logger)

	return data


