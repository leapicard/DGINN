#import DataObject, AnalysisFunc
import logging, os, shlex, sys
from collections import defaultdict, OrderedDict
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SearchIO
from time import sleep

def parseBlast(blastRes):

	with open(blastRes, "r") as blast:
		listBlastRes = blast.readlines()
		blast.close()		
		listAcc = []
		
		for hit in listBlastRes:
			if not hit.startswith("#"):
				listAcc.append(hit.split("|")[-2])
		
	return(listAcc)

def blast(queryFile, output_file, db, evalue, percId, cov, apiKey, remote, query):
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
	blastRes = output_file
	#logger = logging.getLogger("main.blast")
	#logger.info("Running Blast")

	if remote:
		
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
											  hitlist_size=1000, 
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



if __name__ == "__main__" :

    queryFile = snakemake.input[0]
    o = snakemake.output[0]

    #queryFile = "DGINN/examples/ex_CCDS.fasta"
    #o = "blastres.tsv"
    db = "nr"
    evalue = 0.0001
    percId = 70
    cov = 50
    apiKey = ""
    remote = True
    query = ""
    blast(queryFile, 
			      o,
			      db, 
			      evalue, 
			      percId, 
			      cov, 
			      apiKey, 
			      remote, 
			      query)