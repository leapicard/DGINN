#import DataObject, AnalysisFunc
import logging, os, shlex, sys
from collections import defaultdict, OrderedDict
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SearchIO, SeqIO
from time import sleep
import pickle
from AnalysisFunc import cmd

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
	logger = logging.getLogger("main.blast")
	logger.info("Running Blast")

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
		cmd(cmdBlast, False)
		logger.debug("Blast command: {:s}".format(cmdBlast))

	if os.path.exists(blastRes) and os.stat(blastRes).st_size > 0:								
		logger.info("Blast results written to: {:s}, {} sequences retrieved.".format(blastRes, nbSeq))
	else:
		logger.error("Blast didn't run, exiting DGINN.")
		sys.exit()
	
	return(blastRes)

def setGenAttr(data,params):	
		#Set attributes from the queryFile.
		step = params["step"]
		queryFile = data["queryFile"]

		if step == "blast":
			accns = list(SeqIO.parse(open(queryFile),'fasta'))
			accn = accns[0]
			data["queryName"] = accn.id
			params["queryName"] = accn.id
		
		return data,params

if __name__ == "__main__" :
	with open(sys.argv[1], 'rb') as fichier:
		params = pickle.load(fichier)
	with open(sys.argv[2], 'rb') as fichier:
		data = pickle.load(fichier)
	
	o = "results/blastres.tsv"
	
	data["blastRes"] = blast(data["queryFile"], 
			      o,
			      data["db"], 
			      params["evalue"], 
			      params["percID"], 
			      params["mincov"], 
			      params["APIKey"], 
			      params["remote"], 
			      params["entryQuery"])

	data["lBlastRes"] = parseBlast(data["blastRes"])
	
	data,params = setGenAttr(data,params)

	with open(sys.argv[3],'wb') as fichier_data:
		pickle.dump(data,fichier_data,pickle.HIGHEST_PROTOCOL)
	with open(sys.argv[4],'wb') as fichier_data:
		pickle.dump(params,fichier_data,pickle.HIGHEST_PROTOCOL)
	