
#coding: utf8
from Bio import SeqIO
from statistics import median, mean
import fastaResFun

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
	#logger = logging.getLogger("main.alignment")
	
	if queryName in dId2Seq:
		print(f"Discarding sequences with less than {cov} coverage of query.\n")
		print("Discarding sequences with less than {:d}% coverage of query.".format(cov))
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
			outC.write(fastaResFun.dict2fasta(dKeep))
			outC.close()
		print(f"Discarded {nbOut} sequences\n")
		print("Discarded {:d} sequences".format(nbOut))
	
		return(outCov, nbOut)
	
	else:
		#logger.warning("Provided query name not found in the alignment, skipping coverage check.")
		return(aln, 0)