from Bio import SeqIO
from ete3 import Tree

"""The file pools all the functions needed to check formats of the user-provided files at each step of the pipeline."""

def isCCDSFasta(targetFile):
	"""
	Check if provided file is a fasta file containing one CCDS.
	
	@param targetFile: file to be tested
	@return boolFile: boolean of answer
	"""
	accns = list(SeqIO.parse(targetFile,'fasta'))
	boolFile = any(accns)
	seq = ""
	
	for rec in accns:
		seq = str(rec.seq)
		break
	
	if len(accns) > 1:
		boolFile = False
	
	if '-' in seq and not seq.startswith("ATG"):
		boolFile = False
	
	return(boolFile)

def isAccns(targetFile):
	"""
	Check if provided file is at least not a fasta file for accessing the accession step.
	
	@param targetFile: file to be tested
	@return boolFile: boolean of answer
	"""
	boolFile = True
	with open(targetFile, "r") as target:
		bad = any(">" in line for line in target.readlines())
		target.close()		
		if bad:
			boolFile = False
	
	return(boolFile)

def isAln(targetFile):
	"""
	Check if provided file is a fasta file of an alignment (check if all sequences are the same length).
	
	@param targetFile: file to be tested
	@return boolFile: boolean of answer
	"""
	accns = list(SeqIO.parse(targetFile,'fasta'))
	boolFile = any(accns)
	lSeq = [len(rec.seq) for rec in accns]
	
	if not all(l == lSeq[0] for l in lSeq):
		boolFile = False
	
	return(boolFile)	
	
def isTree(targetFile):
	"""
	Check if provided file is a file of a phylogenetic tree.
	
	@param targetFile: file to be tested
	@return boolFile: boolean of answer
	"""
	try:
		boolFile = any(Tree(targetFile))
	except:
		boolFile = False
	
	return(boolFile)

def isBlastRes(targetFile):
	"""
	Check if provided file is a file of a phylogenetic tree.
	
	@param targetFile: file to be tested
	@return boolFile: boolean of answer
	"""
	boolFile = False
	with open(targetFile, "r") as target:
		if target.readline().startswith("# BLAST"):
			boolFile = True
		target.close()	
	return(boolFile)
	
def isFasta(targetFile):
	"""
	Check if provided file is a simple fasta and contain multiple sequences.
	
	@param targetFile: file to be tested
	@return boolFile: boolean of answer
	"""
	
	accns = list(SeqIO.parse(targetFile,'fasta'))
	boolFile = any(accns)
	
	if len(accns) < 2:
		boolFile = False
	
	return(boolFile)	

def isAllORFs(targetFile):
	"""
	Check if provided file is a file of unaligned coding sequences.
	
	@param targetFile: file to be tested
	@return boolFile: boolean of answer
	"""
	accns = list(SeqIO.parse(targetFile,'fasta'))
	boolFile = any(accns)
	lSeq = [rec.seq for rec in accns]
	
	if not all(l.startswith("ATG") for l in lSeq):
		boolFile = False
	
	return(boolFile)
	
def warningCheck(targetFile):
	"""
	Check if enough sequences for proper study (at least 8).
	
	@param targetFile: file to be tested
	@return boolFile: boolean of answer
	"""
	accns = list(SeqIO.parse(targetFile,'fasta'))
	boolFile = any(accns)
	if len(accns) < 8:
		boolFile = False
	return(boolFile)

def isFastaforBlast(targetFile):
	file = list(SeqIO.parse(targetFile,'fasta'))
	if len(file) == 1:
		return True


