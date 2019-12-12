import logging, os, sys
from Bio import SeqIO

"""
This file pools all the object and related functions used in the G2P pipeline.
"""

class basicData:
	"""
	Object pooling all the necessary information for running the pipeline as well as the list of CCDS accessions and Blast result file.
	
	@attr1 QueryFile: Path to the input file
	@attr2 o: Path to the output directory (user determined or default)
	@attr3 logger: Logger object
	@attr4 db: Name of the Blast database
	@attr5 geneName: Gene name
	@attr6 sequence: Gene sequence
	@attr8 blastRes: Path to Blast result file
	@attr9 lBlastRes: homologues of the gene // STILL USEFUL? CHECK
	@attr10 sptree: Path to the species tree
	@attr11 aln: Path to the alignment file
	@attr12 tree: Path to the phylogenetic tree file 
	@attr13 ctrlFile: Path to the Control file for positive selection analyses
	@attr14 alnFormat: Format of alignment files
	@attr15 baseName: base used for naming convention
	"""
	def __init__(self, inFile, outDir, database, timeStamp, spTree, aln, tree, queryName):
		self.queryFile = inFile
		self.o = outDir
		self.logger = logging.getLogger("main")
		self.db = database
		self.geneName = ""
		self.sequence = ""
		self.blastRes = ""
		self.lBlastRes = []
		self.sptree = spTree
		self.ORFs = ""
		self.aln = aln
		self.tree = tree
		self.queryName = queryName
		self.alnFormat = "Fasta"
		self.baseName = ""
		
		## init self.o
		if os.path.exists(self.o):
			self.o = os.path.abspath(self.o)+"/"
		elif self.o != "" and not os.path.exists(self.o):
			try:
				os.makedirs(self.o)
			except FileExistsError:
				pass
			if self.o.endswith("/"):
				self.o = self.o+"/"
		else:
			wrongDir = self.o
			if inFile != "":
				self.o = inFile.split(".")[0]+"_results_"+timeStamp+"/"
			else:
				self.o = aln.split(".")[0]+"_results_"+timeStamp+"/"
			if wrongDir == "":
				try:
					os.makedirs(self.o)
				except FileExistsError:
					pass
			else:
				self.logger.warn("Provided path to output folder {:s} does not exist, defaulting to {:s}".format(wrongDir, self.o))
	
	
	def setGenAttr(self, step):
		
		#Set attributes from the queryFile.
		
		if step == "blast":
			accns = list(SeqIO.parse(open(self.queryFile),'fasta'))
			accn = accns[0]
			self.queryName = accn.id
		else:
			if self.queryName == "":
				self.logger.warn("queryName is not provided, coverage check will not run.")
	
