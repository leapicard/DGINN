import G2P_Object, AnalysisFunc, LoadFileFunc
import logging, subprocess, shlex
from Bio import SeqIO
from multiprocessing import Pool
from collections import defaultdict, OrderedDict
from statistics import median
import requests
from ete3 import NCBITaxa

"""
This file pools functions related to the creation and conversion of fasta format data.
"""

def dict2fasta(dico):
	"""
	Function converting a dictionary associating gene names to their sequences into a writable Fasta format.

	@param dico: Dictionary associating gene names (keys) to their CCDS fasta sequence (values)
	@return txtoutput: Fasta formatted string of the dictionary
	"""

	txtoutput = ""
	for key, value in dico.items():
		txtoutput += ">{:s}\n{:s}\n".format(str(key),str(value))

	return(txtoutput)

def blastExtract(accnFile, db, lBlastRes, geneName, geneDir, apiKey, remote):
	"""
	Function downloading sequences for a list of genes using their accessions.

	@param1 accnFile: Path
	@param2 db: Database name
	@param3 lBlastRes: List of accessions
	@param4 geneName: Gene name
	@param5 geneDir: Path directory 
	@param6 apiKey: Key for the API of NCBI
	@param7 remote: Boolean (online database or not)
	@return out: Path to the file containing the extracted sequences
	"""

	## extract fasta from blast database

	out = geneDir+accnFile.replace("_accns.txt", "_primatesnoquery.fasta").split("/")[-1]
	logger = logging.getLogger("main")

	#If remote is set to False, we use the database in local
	if remote == "False":
		AnalysisFunc.cmd("blastdbcmd -db {:s} -dbtype nucl -entry_batch {:s} -out {:s} -outfmt %f".format(db, accnFile, out), False)
	#If remote is set to True and the lengh of lBlastRes isn't big, generate one link
	elif len(lBlastRes) < 20:
		if apiKey != "":
			link = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sequences&id={:s}&rettype=fasta&api_key={:s}".format(",".join(lBlastRes), apiKey)
		else:
			link = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sequences&id={:s}&rettype=fasta".format(",".join(lBlastRes))
		r = requests.get(link, headers={"Content-Type" : "text/x-fasta"})
		fasta = open(out, "w")
		texTemp = r.text.split("\n")
		text = [ i for i in texTemp if i != "" ]
		fasta.write("\n".join(text))
		fasta.close()
	#If there're too much accessions number, generate several link to get fasta
	else:
		start = 0
		end = 20
		lRes = []
		while end < len(lBlastRes):
			if apiKey != "":
				link = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sequences&id={:s}&rettype=fasta&api_key={:s}".format(",".join(lBlastRes[start:end]), apiKey)
			else:
				link = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sequences&id={:s}&rettype=fasta".format(",".join(lBlastRes[start:end]))
			r = requests.get(link, headers={"Content-Type" : "text/x-fasta"})
			texTemp = r.text.split("\n")
			text = [ i for i in texTemp if i != "" ]
			lRes.append("\n".join(text))
			start=end
			end +=20
		if len(lBlastRes) != end:
			if apiKey != "":
				link = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sequences&id={:s}&rettype=fasta&api_key={:s}".format(",".join(lBlastRes[start:]), apiKey)
			else:
				link = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sequences&id={:s}&rettype=fasta".format(",".join(lBlastRes[start:]))
			r = requests.get(link, headers={"Content-Type" : "text/x-fasta"})
			texTemp = r.text.split("\n")
			text = [ i for i in texTemp if i != "" ]
			lRes.append("\n".join(text))
			
		#Filtre à ajouter

		with open(out, "w") as fasta:
			fasta.write("\n".join(lRes))
		
	#logger.debug(cmd)
	logger.info("{:s}: extracted fasta sequences from blast results.".format(geneName.split("|")[1]))
	
	return(out)

def catFile(lBlastRes, geneName, sequence, hitsFasta, sptree, o, apiKey, treerecs, logger):
	"""
	Function dowloading species to generate new datas.

	@param1 lBlastRes: List of accessions
	@param2 geneName: Gene name
	@param3 sequence: gene sequence
	@param4 hitsFasta: Path
	@param5 sptree: species tree
	@param6 o: Output directory
	@param7 apiKey: Key for the API of NCBI
	@param8 treerecs: Boolean
	@param9 logger: Object logging
	@return1 outCat: Path to the file containing the sequences and the new IDs
	@return2 corSG: Path
	"""
	dSpecies = {}
	for i in lBlastRes:

		if apiKey != "":
			link = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sequences&id={:s}&api_key={:s}&retmode=text".format(i, apiKey)
		else:
			link = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sequences&id={:s}&retmode=text".format(i)
		r = requests.get(link, headers={"Content-Type" : "text"})
		texTemp = r.text.split("\n")

		line = 0
		while "taxname" not in texTemp[line]:
			line += 1
		dSpecies[i] = texTemp[line].split('"')[1]

	dNewId2Seq = OrderedDict()
	dNewId2Len = {}

	lName = geneName.split("|")
	newID = lName[0]+lName[2]
	dNewId2Seq[newID] = sequence
	for accn in SeqIO.parse(open(hitsFasta),'fasta'):
		if accn.id in dSpecies:
			accnNb, accnSeq = accn.id, str(accn.seq)
			spInter = dSpecies[accn.id].split(" ")
			accnSp = spInter[0][:3].lower()+"".join([ i[:3].title() for i in spInter[1:]])

			if "." in accnNb:
				accnNb = accnNb.replace(".", "dot")

			accnNewID = accnSp+accnNb
			dNewId2Seq[accnNewID] = accnSeq
			dNewId2Len[accnNewID] = len(accnSeq)

	logger.debug(dNewId2Seq)
	logger.debug(dNewId2Len)		

	m = median(dNewId2Len.values())
	for k, v in dNewId2Len.items():
		if v > 10000 and v > 5*m:
			try:
				del dNewId2Seq[k]
				logger.debug("{:s}: deleted sequence {:s} (length {:d})".format(geneName, k, v))
			except KeyError:
				pass
	
	outCat = hitsFasta.replace("noquery", "")
	with open(outCat, "w") as out:
		out.write(dict2fasta(dNewId2Seq))
	
	logger.info("{:s}: added human sequence to file and deleted sequences with excessive length.".format(geneName))
	
	corSG=""
	if treerecs == "True":
		outCat, corSG = LoadFileFunc.filtreData(sptree, outCat, o)

	return(outCat, corSG)


def fastaCreation(data, logger, remote, apiKey, treerecs):
	"""
	Function handling the creation of fasta files in the pipeline.

	@param1 data: basicdata object
	@param2 logger: Logger object
	@param3 remote: Boolean (online database or not)
	@param4 apiKey: Key for the API of NCBI
	@param5 treerecs: Booleans
	"""
	listExtract = blastExtract(data.accnFile, data.db, data.lBlastRes, data.geneName, data.geneDir, apiKey, remote)
	setattr(data, "hitsFasta", listExtract)
	logger.info("Extracted blast hits from blast database")

	listCat, corSG = catFile(data.lBlastRes, data.geneName, data.sequence, data.hitsFasta, data.sptree, data.o, apiKey, treerecs, logger)
	setattr(data, "catFile", listCat)
	setattr(data, "cor", corSG)
	logger.info("Added human sequence to fasta file")
