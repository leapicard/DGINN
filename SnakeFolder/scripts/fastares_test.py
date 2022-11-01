import os
import loadfile_test, tree_test
from blast_test import parseBlast
import logging, sys
from Bio import Entrez
from statistics import median

def dict2fasta(dico):

	txtoutput = ""
	for key, value in dico.items():
		txtoutput += ">{:s}\n{:s}\n".format(str(key),str(value))

	return(txtoutput)


def remoteDl(lBlastRes, queryName, apiKey):
	dId2Seq = {}
	lTax = []

	Entrez.email="example@example.com"
	Entrez.api_key = apiKey
	
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
                                name = "pot"+queryName.split("_")[1]
			except IndexError:
                                name = "pot"
		if tax == "synCon" or 'GBSeq_sequence' not in record.keys():
			continue
		else:
			lTax.append(tax)
			dId2Seq[tax+"_"+name+"_"+acc.split(".")[0]] = record['GBSeq_sequence'].upper()
			
	handle.close()
	
	return(dId2Seq)


def sizeCheck(dId2Seq):
	
	dId2Len = {Id:len(seq) for Id, seq in dId2Seq.items()}
	m = median(dId2Len.values())
	n = 0
	
	for k, v in dId2Len.items():
		if m > 10000:
			if v > 2*m:
				try:
					del dId2Seq[k]
					n += 1
				except KeyError:
					pass
		else:
			if v > 3*m or v > 20000:
				try:
					del dId2Seq[k]
					n += 1
				except KeyError:
					pass
	
	return(dId2Seq)
	

def catFile(queryFile, dId2Seq, firstFasta):
	
	with open(queryFile, "r") as query:
		lquery = query.readlines()
		query.close()
		dId2Seq[lquery[0].strip().replace(">", "")] = lquery[1]
	
		with open(firstFasta, "w") as fasta:
			fasta.write(dict2fasta(dId2Seq))
			fasta.close()	
	return(firstFasta)
	

def fastaCreation(queryFile, tsvfile, remote, apiKey, maxLen, step, query, treerecs):
	lBlastRes = parseBlast(tsvfile)

	if remote:
		dId2Seq = remoteDl(lBlastRes,"", apiKey)
	
	dId2Seq = sizeCheck(dId2Seq)
	
	firstFasta = snakemake.output[0]

	if step == "blast":
		firstFasta = catFile(queryFile, dId2Seq, firstFasta)
	else:
		with open(firstFasta, "w") as out:
			out.write(dict2fasta(dId2Seq))
			out.close()

	if treerecs:
		outCat, corSG = loadfile_test.filterData(snakemake.input[2], firstFasta, "results/")


if __name__ == "__main__" :
	tsvfile = snakemake.input[0]
	queryFile = snakemake.input[1]
	remote = True
	apiKey = ""
	maxLen = ""
	step = ""
	query = ""
	treerecs = True
	fastaCreation(queryFile,tsvfile, remote, apiKey, maxLen, step, query, treerecs)
