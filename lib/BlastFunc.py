import DataObject, AnalysisFunc
import logging, os, shlex, sys
from collections import defaultdict, OrderedDict
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SearchIO
from time import sleep

"""This file pools the necessary functions to run Blast and parse its results."""

def blast(queryFile, outDir, baseName, db, evalue, percId, cov, apiKey, remote, query):
	"""
	Function running Blast on the provided gene list.
	
	@param1 queryFile: Path
	@param2 db: Database name
	@param3 evalue: Blast e-value
	@param4 percId: Blast percentage of identity
	@param5 cov: Blast coverage
	@param6 apiKey: Key for the API of NCBI
	@param7 remote: Boolean (online database or not)
	@param8 query: Strings which will select the database online
	@return blastRes: Path to Blast results file
	"""
	# Blast results file
	blastRes = outDir+baseName+"_blastres.tsv"
	logger = logging.getLogger("main.blast")
	logger.info("Running Blast")

	if remote:
		
		if query == "":
			logger.info("Database {:s} will be searched for all species".format(db))
		else:
			logger.info("Database search on {:s} will be limited to {:s}".format(db, query))
			
		sequence = open(queryFile).read()
		seqL = len(sequence.split("\n")[1])
		url = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi?api_key={:s}'.format(apiKey)
		
		blasted = False
		while not blasted:
			try:
				sleep(2)
				resultHandle = NCBIWWW.qblast("blastn", 
											  db, 
											  sequence, 
											  entrez_query = query, 
											  url_base = url, 
											  hitlist_size=100000, 
											  perc_ident = percId, 
											  threshold = evalue)
				blasted = True
			except ValueError:
				blasted = False
			
		f = open(blastRes, "w")
		nbSeq = 0

		for record in NCBIXML.parse(resultHandle):
			#f.write("# "+record.application+" "+record.version+"\n# Query: "+record.query+"\n# Database: "+record.database+"\n# Fields: subject id, ?\n# "+str(len(record.alignments))+" hits found\n")
			f.write("# {} {}\n# Query: {}\n# Database: {}\
					\n# Fields: subject id, ?\n# {} hits found\n".format(record.application,
																		 record.version,
																		 record.query,
																		 record.database,
																		 str(len(record.alignments))))
			
			for alignment in record.alignments:
				for hsp in alignment.hsps:
					qcov = hsp.align_length/seqL*100
					if qcov > cov:
						f.write(alignment.title+"\n")
						nbSeq += 1
		f.close()

	else:
		cmdBlast = "blastn -task blastn -db {:s} -query {:s} -out {:s} \
					-evalue {:f} -perc_identity {:d} -qcov_hsp_perc {:d} \
					-outfmt \"7 sseqid\" -max_target_seqs 1000".format(db, 
																	   queryFile, 
																	   blastRes, 
																	   evalue, 
																	   percId, 
																	   cov)
		AnalysisFunc.cmd(cmdBlast, False)
		logger.debug("Blast command: {:s}".format(cmdBlast))

	if os.path.exists(blastRes) and os.stat(blastRes).st_size > 0:								
		logger.info("Blast results written to: {:s}, {} sequences retrieved.".format(blastRes, nbSeq))
	else:
		logger.error("Blast didn't run, exiting DGINN.")
		sys.exit()

	return(blastRes)

def parseBlast(blastRes):
	"""
	Function parsing accessions of Blast results per gene.

	@param1 blastRes: Path
	@param2 logger: Logging object
	@return lGeneRes: list of all the accessions of homologous genes found through Blast
	"""
	logger = logging.getLogger("main.accessions")
	
	with open(blastRes, "r") as blast:
		listBlastRes = blast.readlines()
		blast.close()		
		listAcc = []
		
		for hit in listBlastRes:
			if not hit.startswith("#"):
				listAcc.append(hit.split("|")[-2])
		
	logger.debug(listAcc)
	logger.info("Extracted gene accessions from blast results")
	return(listAcc)

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
	@return data: List of gene objects
	"""

	data.blastRes = blast(data.queryFile, 
			      data.o, 
			      data.baseName, 
			      data.db, 
			      evalue, 
			      percId, 
			      cov, 
			      apiKey, 
			      remote, 
			      query)
	data.lBlastRes = parseBlast(data.blastRes)

	return data


