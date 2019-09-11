import G2P_Object, AnalysisFunc, LoadFileFunc
import logging, subprocess, shlex, sys, requests, json
from Bio import SeqIO, Entrez
from collections import defaultdict, OrderedDict
from statistics import median
from ete3 import NCBITaxa
from time import sleep

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

def blastExtract(accnFile, db, lBlastRes, geneName, outDir, apiKey, remote):
	"""
	Function downloading sequences for a list of genes using their accessions.

	@param1 accnFile: Path
	@param2 db: Database name
	@param3 lBlastRes: List of accessions
	@param4 geneName: Gene name
	@param5 outDir: Path directory 
	@param6 apiKey: Key for the API of NCBI
	@param7 remote: Boolean (online database or not)
	@return out: Path to the file containing the extracted sequences
	"""

	## extract fasta from blast database

	out = outDir+accnFile.replace("_accns.txt", "_primatesnoquery.fasta").split("/")[-1]
	logger = logging.getLogger("main")

	#If remote is set to False, we use the database in local
	if not remote:
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
		
	logger.info("Extracted fasta sequences from blast results: {:s}".format(out))
	
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
	
	"""
	Entrez.api_key = apiKey
	Entrez.email = ""
	"""
	for i in lBlastRes:
		
		if apiKey != "":
			link = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sequences&id={:s}&api_key={:s}&retmode=text".format(i, apiKey)
		else:
			link = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sequences&id={:s}&retmode=text".format(i)
	
		r = requests.get(link)
		while r.status_code != 200:
			sleep(2)
			r = requests.get(link)
		
		handle = r.text
		"""	
		handle = Entrez.efetch(db = "nucleotide", id = i, retmode = "text")
		handle = handle.read()
		"""
		tax = handle.split('taxname')[1].split('"')[1].split(" ")
		tax = tax[0][:3].lower()+"".join([ i[:3].title() for i in tax[1:]])
		

		try:
			name = handle.split('locus')[1].split('"')[1].upper().replace('/', "_").replace(" ", "").replace("-", "")
			#print(name)
			if "{" in name:
				name = "pot"+geneName.split("|")[1]
			if " " in name:
				name = "pot"+geneName.split("|")[1]
			if "\n" in name:
				name = "pot"+geneName.split("|")[1]
		except:
			name = "pot"+geneName.split("|")[1]
			
		dSpecies[i] = tax+"_"+name
	
	dNewId2Seq = OrderedDict()
	dNewId2Len = {}

	newID = geneName.replace("|", "_")
	dNewId2Seq[newID] = sequence
	for accn in SeqIO.parse(open(hitsFasta),'fasta'):
		if accn.id in dSpecies:
			accnNb, accnSeq = accn.id, str(accn.seq)
			accnSp = dSpecies[accn.id]
			"""spInter = dSpecies[accn.id].split(" ")
			accnSp = spInter[0][:3].lower()+"".join([ i[:3].title() for i in spInter[1:]])"""

			if "." in accnNb:
				accnNb = accnNb.replace(".", "dot")

			accnNewID = accnSp+"_"+accnNb
			dNewId2Seq[accnNewID] = accnSeq
			dNewId2Len[accnNewID] = len(accnSeq)

	#logger.debug(dNewId2Seq)
	#logger.debug(dNewId2Len)		

	m = median(dNewId2Len.values())
	for k, v in dNewId2Len.items():
		if m > 10000:
			if v > 2*m:
				try:
					del dNewId2Seq[k]
					logger.debug("Deleted sequence {:s} (length {:d})".format(k, v))
				except KeyError:
					pass
		else:
			if v > 3*m or v > 20000:
				try:
					del dNewId2Seq[k]
					logger.debug("Deleted sequence {:s} (length {:d})".format(k, v))
				except KeyError:
					pass
	
	outCat = hitsFasta.replace("noquery", "")
	with open(outCat, "w") as out:
		out.write(dict2fasta(dNewId2Seq))
	
	logger.info("Added original query sequence to file and deleted sequences with excessive length: {:s}".format(outCat))
	
	return(outCat)


def fastaCreation(data, logger, remote, apiKey, treerecs):
	"""
	Function handling the creation of fasta files in the pipeline.

	@param1 data: basicdata object
	@param2 logger: Logger object
	@param3 remote: Boolean (online database or not)
	@param4 apiKey: Key for the API of NCBI
	@param5 treerecs: Booleans
	"""
	listExtract = blastExtract(data.accnFile, data.db, data.lBlastRes, data.geneName, data.o, apiKey, remote)
	setattr(data, "hitsFasta", listExtract)

	listCat = catFile(data.lBlastRes, data.geneName, data.sequence, data.hitsFasta, data.sptree, data.o, apiKey, treerecs, logger)
	if treerecs:
		outCat, corSG = LoadFileFunc.filterData(data.sptree, listCat, data.o)
		setattr(data, "catFile", outCat)
		setattr(data, "cor", corSG)
