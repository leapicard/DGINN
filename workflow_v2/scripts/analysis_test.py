import sys
import logging, subprocess, shlex, os, ete3
from Bio import SeqIO, AlignIO
from collections import defaultdict, OrderedDict
from itertools import chain
from statistics import median, mean
import pandas as pd
import pickle
import fastares_test

def cmd(commandLine, choice, verbose = False):
	"""
	Function executing a command line in a bash terminal.

	@param1 commandLine: String corresponding to a bash command line
	@param2 choice: Boolean determining whether the command is executed within the shell 
	"""
    if verbose :
        stdout=None
    else :
        stdout=subprocess.PIPE
    
    lCmd = shlex.split(commandLine)
    try:
	    run = subprocess.call(lCmd, 
			        shell=choice,
                          stdout=stdout,
			        stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as err:
        sys.stderr.write(str(err))


def getORFs(catFile, queryName, geneDir):
	"""
	Function to find Open Reading Frames within the sequence of each gene and select the longest one.

	@param1 catFile: Path
	@param2 geneName: Gene name
	@param3 geneDir: Gene directory
	@return outORF: Path to the file containing the longest ORFs
	"""

	outORFraw = geneDir
	logger = logging.getLogger("main.orf")
	
	logger.debug("getorf -sequence {:s} -outseq {:s} -table 0 -find 3 -noreverse".format(catFile, outORFraw))
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
			logger.debug("Deleted sequence {:s} (duplicate)".format(i))
		
	logger.info("Deleted {} sequences as duplicates".format(n))
	
	outORF = sys.argv[4]
	
	a = fastares_test.dict2fasta(dId2Longest)

	with open(outORF, "w") as outO:
	    outO.write(fastares_test.dict2fasta(dId2Longest))
        outO.close()
	  
    logger.info("Extracted longest ORFs: {:s}".format(outORF))

    return(outORF)


def orfFinder(data):
    """
	Procedure which launch the ORF step

	@param1 data: basicdata object
	@param2 logger: Logging object
	"""

    ORFile = getORFs(data["seqFile"],
					 data["queryName"], 
					 sys.argv[2])
    data["ORFs"] = ORFile


    

if __name__ == "__main__" :
    with open(sys.argv[1], 'rb') as fichier:
        data = pickle.load(fichier)

    data = orfFinder(data)
    
    with open(sys.argv[3],'wb') as fichier_data:
        pickle.dump(data,fichier_data,pickle.HIGHEST_PROTOCOL)
    