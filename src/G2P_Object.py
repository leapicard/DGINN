import logging, os
from Bio import SeqIO

"""
This file pools all the object and related functions used in the G2P pipeline.
"""

class basicData:
	"""
	Object pooling all the necessary information for running the pipeline as well as the list of CCDS accessions and Blast result file.
	
	@attr1 CCDSFile: Path to the input file
	@attr2 o: Path to the output directory (user determined or default)
	@attr3 logger: Logger object
	@attr4 db: Name of the Blast database
	@attr5 geneName: Name's gene
	@attr6 sequence: sequence's gene
	@attr8 blastRes: Path to Blast result file for the provided gene list
	@attr9 lBlastRes: homolgues of the gene
	@attr10 sptree: Path of the species tree
	@attr11 aln: Path of the alignment file
	@attr12 tree: Path of the tree file 
	@attr13 ctrlFile: Path of the Control file for PSP step
	@attr14 alnFormat: Format of alignments files
	@attr15 baseName: Base to name files
	"""
	def __init__(self, inFile, outDir, database, timeStamp, spTree, aln, tree):
		self.CCDSFile= inFile
		self.o = outDir
		self.logger = logging.getLogger("main")
		self.db = database
		self.geneName = ""
		self.sequence = ""
		self.blastRes = ""
		self.lBlastRes = []
		self.sptree = spTree
		self.aln = aln
		self.tree = tree
		self.alnFormat = "Fasta"
		self.baseName = ""
		## initiate self.o
		if os.path.exists(self.o):
			self.o = os.path.abspath(self.o)+"/"
		else:
			wrongDir = self.o
			if inFile != "":
				self.o = inFile.split(".")[0]+"_results_"+timeStamp+"/"
			else:
				self.o = aln.split(".")[0]+"_results_"+timeStamp+"/"
			if wrongDir == "":
				os.makedirs(self.o)
			else:
				self.logger.warn("Provided path to output folder {:s} does not exist, defaulting to {:s}".format(wrongDir, self.o))

	def setGenAttr(self):
		"""
		Set attributs from the CCDSFile.
		"""
		accns = list(SeqIO.parse(open(self.CCDSFile),'fasta'))
		accn = accns[0]
		self.geneName = accn.id
		self.sequence = str(accn.seq)
		if len(accns) > 1:
			self.logger.info("Too many sequences; only the first one will be treated...")
