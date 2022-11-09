
#coding: utf8
import sys
import logging, subprocess, shlex, os, ete3
from Bio import SeqIO, AlignIO
from collections import defaultdict, OrderedDict
from itertools import chain
from statistics import median, mean
import pandas as pd

import pickle
import fastares_test

def covAln(aln, cov, queryName, o):
	"""
	Function to discard sequences from alignment according to coverage to query.

	@param1 aln: Path to prank alignment
	@param2 cov: minimum coverage necessary to keep sequence (from parameter file)
	@param4 queryName: full identifier of the sequence against which to check coverage
	@param5 o: output directory
	@return outCov: Path to file of conserved sequences
	"""
	
	dId2Seq = {fasta.id:str(fasta.seq) for fasta in SeqIO.parse(open(aln),'fasta')}
	logger = logging.getLogger("main.alignment")
	
	if queryName in dId2Seq:
		logger.info("Discarding sequences with less than {:d}% coverage of query.".format(cov))
		outCov = o
		
		nbOut = 0
		lIndexes = [pos for pos, char in enumerate(dId2Seq[queryName]) if char != "-"]
		

		dKeep = {}
		for ID, seq in dId2Seq.items():
			seqPos = [seq[x] for x in lIndexes]
			seqCov = (len(seqPos) - seqPos.count("-"))/len(seqPos)*100
			
			if seqCov > cov:
				dKeep[ID] = seq
		
		nbOut = len(dId2Seq) - len(dKeep)
		
		with open(outCov, "w") as outC:
			outC.write(fastares_test.dict2fasta(dKeep))
			outC.close()
		logger.info("Discarded {:d} sequences".format(nbOut))
	
		return(outCov, nbOut)
	
	else:
		logger.warning("Provided query name not found in the alignment, skipping coverage check.")
		return(aln, 0)




if __name__ == "__main__" :

	"""
	arg1 = mafft output : fasta alignment file path
	arg2 = covAln output : fasta alignment file path
	arg3 = data dictionnary : pickle file containing dictionnary of dataObject path, before mafft update
	arg4 = data dictionnary : pickle file containing dictionnary of dataObject path, after mafft update
	"""


	with open(sys.argv[3], 'rb') as fichier:
		data = pickle.load(fichier)

	data["aln"] = sys.argv[1]
	fasCov, nbOut = covAln(aln = data["aln"], 
							      cov = 50, 
							      queryName = data["queryName"], 
							      o = sys.argv[2])
    
	with open(sys.argv[4],'wb') as fichier_data:
		pickle.dump(data,fichier_data,pickle.HIGHEST_PROTOCOL) 