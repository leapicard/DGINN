import LoadFileFunc, TreeFunc
import logging
import sys
from Bio import Entrez
from statistics import median

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
        txtoutput += ">{:s}\n{:s}\n".format(str(key), str(value))

    return txtoutput



def remoteDl(lBlastRes, queryName, apiKey):
	"""
	Function dowloading species to generate new datas.

	@param1 lBlastRes: List of accessions
	@param2 queryName
	@param3 apiKey: apikey, not considered if equals ""
	@return dId2Seq: dictionnary from gene ids to sequences
	"""
	logger = logging.getLogger("main.accessions")
	dSpecies = {}
	dId2Seq = {}
	lTax = []

	Entrez.email="example@example.com"
	if apiKey!="":
	  Entrez.api_key = apiKey
	  logger.info("No ApiKey")
	else:
	  logger.info("ApiKey " + apiKey)
          
	handle = Entrez.efetch(db="nuccore", id=lBlastRes , idtype="acc", retmode="xml")
	records = list(Entrez.read(handle))
	
	for record in records:
		acc = record['GBSeq_primary-accession']
		tax = record['GBSeq_organism']
		if " x " in tax:
			tax = tax.split(" x ")[0].split(" ")
		elif " X " in tax:
			tax = tax.split(" X ")[0].split(" ")
		else:
			tax = tax.split(" ")
			
		tax = tax[0][:3].lower()+"".join([ i[:3].title() for i in tax[1:]])

		features = [record['GBSeq_feature-table'][i]['GBFeature_quals'] for i, d in enumerate(record['GBSeq_feature-table']) if 'GBFeature_quals' in d]
		
		for feat in features:
			for d in feat:
				if ('GBQualifier_name', 'gene') in d.items():
					name = d['GBQualifier_value']
					break
				else:
					name = ""
		if "." in name or "-" in name or name == "":
			try:
				name = "pot_"+queryName.split("_")[1]
			except IndexError:
				name = "pot"
		if tax == "synCon" or 'GBSeq_sequence' not in record.keys():
			continue
		else:
			lTax.append(tax)
			dId2Seq[tax+"_"+name+"_"+acc.split(".")[0]] = record['GBSeq_sequence'].upper()
			
	handle.close()
	nbSp = len(set(lTax))
	print(f"Remote option on, downloaded gene IDs and sequences from NCBI databases ({nbSp} different species represented in the retrieved sequences).\n")
	print("Remote option on, downloaded gene IDs and sequences from NCBI databases ({} different species represented in the retrieved sequences).".format(nbSp))
	
	return(dId2Seq)


def sizeCheck(dId2Seq):
	#logger = logging.getLogger("main.accessions")

	dId2Len = {Id:len(seq) for Id, seq in dId2Seq.items()}
	m = median(dId2Len.values())
	n = 0
	
	for k, v in dId2Len.items():
		if m > 10000:
			if v > 2*m:
				try:
					del dId2Seq[k]
					#logger.debug("Deleted sequence {:s} (length {:d})".format(k, v))
					n += 1
				except KeyError:
					pass
		else:
			if v > 3*m or v > 20000:
				try:
					del dId2Seq[k]
					#logger.debug("Deleted sequence {:s} (length {:d})".format(k, v))
					n += 1
				except KeyError:
					pass
	print(f"Deleted {n} sequences due to excessive length\n")
	print("Deleted {} sequences due to excessive length.".format(n))

	return(dId2Seq)
	

def catFile(queryFile, dId2Seq, firstFasta):
	#logger = logging.getLogger("main")

	with open(queryFile, "r") as query:
		lquery = query.readlines()
		query.close()
		dId2Seq[lquery[0].strip().replace(">", "")] = lquery[1]
	
		with open(firstFasta, "w") as fasta:
			fasta.write(dict2fasta(dId2Seq))
			fasta.close()

	return(firstFasta)
	

def fastaCreation(parameters, lBlastRes, outputfile):
        """
        Function handling the creation of fasta files in the pipeline.
        
        @param1 data: basicdata object
        @param2 remote: Boolean (online database or not)
        @param3 apiKey: Key for the API of NCBI
        @param4 treerecs: Booleans
        """
        
        if parameters["remote"]:
        	dId2Seq = remoteDl(lBlastRes, parameters["queryName"], parameters["APIKey"])
        else: ### need to code this!!!!
        	#logger = logging.getLogger("main.fasta")
        	print("Local retrieval of information not yet implemented, exiting DGINN.")
        	sys.exit()
        
        dId2Seq = sizeCheck(dId2Seq)
        
        firstFasta = outputfile
        
        if parameters["step"] == "blast":
        	firstFasta = catFile(parameters["queryFile"], dId2Seq, firstFasta)
        else:
        	with open(firstFasta, "w") as out:
        		out.write(dict2fasta(dId2Seq))
        		out.close()
        
        if "sptree" in parameters and parameters["duplication"]:
          parameters["sptree"], parameters["duplication"] = TreeFunc.treeCheck(parameters["sptree"], firstFasta)
        
        if parameters["duplication"] and parameters["sptree"]!="":
          parameters["seqFile"], parameters["cor"] = LoadFileFunc.filterData(parameters["sptree"], firstFasta, parameters["outdir"])
	
        return firstFasta
