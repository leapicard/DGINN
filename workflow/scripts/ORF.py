import sys
import logging, subprocess, shlex, os, ete3
from Bio import SeqIO, AlignIO
from collections import defaultdict, OrderedDict
from itertools import chain
from statistics import median, mean
import pandas as pd
import json
import fastaRes, loadFile


def cmd(commandLine, choice, verbose = False):
	"""
	Function executing a command line in a bash terminal.

	@param1 commandLine: String corresponding to a bash command line
	@param2 choice: Boolean determining whether the command is executed within the shell 
	"""
	if verbose:
		stdout=None
	else:
		stdout=subprocess.PIPE
    
	lCmd = shlex.split(commandLine)
	try:
		run = subprocess.call(lCmd, 
			        shell=choice,
                          stdout=stdout,
					stderr=subprocess.PIPE)
	except subprocess.CalledProcessError as err:
	  sys.stderr.write(str(err))


def getORFs(catFile, queryName, outORFraw):
	"""
	Function to find Open Reading Frames within the sequence of each gene and select the longest one.

	@param1 catFile: Path
	@param2 geneName: Gene name
	@param3 outORFraw: path to output file of getORFs.
	@return outORF: Path to the file containing the longest ORFs
	"""

	#logger = logging.getLogger("main.orf")
	
	#logger.debug("getorf -sequence {:s} -outseq {:s} -table 0 -find 3 -noreverse".format(catFile, outORFraw))
	cmd("getorf -sequence {:s} -outseq {:s} -table 0 -find 3 -noreverse".format(catFile, outORFraw), False)
	
	dId2ORFs = defaultdict(list)
	f = SeqIO.parse(open(outORFraw),'fasta')
	for fasta in f:
		fname, fseq = fasta.id, str(fasta.seq)
		if len(fname.split("_")) > 2:
			fname2 = "_".join(fname.split("_")[0:-1])
		else:
			fname2 = fname.split("_")[0]
		dId2ORFs[fname2].append(fseq)
	
	dId2Longest = {}
	for k, v in dId2ORFs.items():
		dId2Longest[k] = max(v, key=len)
		
	# delete duplicate sequences
	dRev = {}
	for k, v in dId2Longest.items():
		dRev.setdefault(v, set()).add(k)
		
	AllDupl = [values for key, values in dRev.items() if len(values) > 1]
	n = 0
	for dupl in AllDupl:
		species = set([x.split("_")[0] for x in dupl])
		
		for sp in species:
			if queryName in dupl:
				firstOcc = queryName
			else:
				lOcc = [x for x in dupl if sp in x]
				
				if len(lOcc) > 0:
					firstOcc = lOcc[0]
				else:
					firstOcc = str(lOcc)
					
			dupl.remove(firstOcc)
		
		for i in dupl:
			dId2Longest.pop(i, None)
			n += 1
			#logger.debug("Deleted sequence {:s} (duplicate)".format(i))
	print(f"Deleted {n} sequences as duplicates\n")
	#logger.info("Deleted {} sequences as duplicates".format(n))
	
	outORF = sys.argv[3]
	#print(fastaRes.dict2fasta(dId2Longest))
	with open(outORF, "w") as outO:
		outO.write(fastaRes.dict2fasta(dId2Longest))
		outO.close()
	
	#logger.info("Extracted longest ORFs: {:s}".format(outORF))

	return(outORF)


def orfFinder(data_dict):
	"""
	Procedure which launch the ORF step

	@param1 config_dict: basicconfig_dict object
	@param2 logger: Logging object
	"""

	ORFile = getORFs(data_dict["seqFile"],
					 data_dict["queryName"], 
					 sys.argv[2])
	data_dict["ORFs"] = ORFile
	return data_dict

    

if __name__ == "__main__" :

	with open(sys.argv[1], 'r') as config_in:
		config_dict = json.load(config_in)

	parameters = config_dict["parameters"]
	data = config_dict["data"]
	data["firstStep"] = "orf" # a enlever après 

	if parameters["step"] == "orf":
		data = loadFile.orfEntry(data)

	if data["firstStep"] == "orf":
		loadFile.spTreeCheck(data, data["firstStep"], parameters["duplication"])

	data = orfFinder(data) 

	config_dict["parameters"] = parameters
	config_dict["data"] = data
    
	with open(sys.argv[1],'w') as config_out:
		json.dump(config_dict, config_out, indent="")
		
		
    