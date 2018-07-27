import G2P_Object, FastaResFunc, AccessionFunc
import logging, subprocess, shlex, os
from Bio import SeqIO, AlignIO
from collections import defaultdict, OrderedDict

######Functions=============================================================================================================

def cmd(commandLine,choice):
	"""
	Function executing a command line in a bash terminal.

	@param1 commandLine: String corresponding to a bash command line
	@param2 choice: Boolean determining whether the command is executed within the shell 
	"""
	
	lCmd = shlex.split(commandLine)
	run = subprocess.run(lCmd, shell=choice, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

######ORF===================================================================================================================

def getORFs(catFile, geneName, geneDir):
	"""
	Function to find Open Reading Frames within the sequence of each gene and select the longest one.

	@param1 catFile: Path
	@param2 geneName: Gene name
	@param3 geneDir: Gene directory
	@return outORF: Path to the file containing the longest ORFs
	"""

	outORFraw = geneDir+catFile.split("/")[-1].split(".")[0]+"_allORFs.fasta"
	print(outORFraw)
	logger = logging.getLogger("main")
	
	cmd("getorf -sequence {:s} -outseq {:s} -table 0 -find 3 -noreverse".format(catFile, outORFraw), False)
	
	#logger.debug(cmd)
	
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
	
	outORF = outORFraw.replace("_allORFs.fasta","_longestORFs.fasta")

	with open(outORF, "w") as outO:
		outO.write(FastaResFunc.dict2fasta(dId2Longest))
	
	logger.info("{:s}: got longest ORFs.".format(geneName.split("|")[1]))

	return(outORF)


def orfFinder(data, logger):
	"""
	Procedure which launch the ORF step

	@param1 data: basicdata object
	@param2 logger: Logging object
	"""
	
	ORFile = getORFs(data.catFile, data.geneName, data.geneDir)
	setattr(data, "ORFs", ORFile)
	logger.info("Got longest ORFs")

#######=================================================================================================================
######PRANK=============================================================================================================
	
def runPrank(ORFs, geneName, outPrank):
	"""
	Function to run PRANK software for codon alignment (Löytynoja, 2014).
	
	@param1 ORFs: Path
	@param2 geneName: Gene name
	@param3 geneDir: Gene Directory
	@return outPrank: Path to Prank results file
	"""
	logger = logging.getLogger("main")
	logger.info("{:s}: start Prank codon alignment.".format(geneName))
	
	cmd("prank -d={:s} -o={:s} -codon -F".format(ORFs, outPrank), False)
	
	logger.info("{:s}: finished Prank codon alignment.".format(geneName.split("|")[1]))


def alnPrank(data, logger):
	"""
	Function creating alignment (aln) attribute in each gene object of the list.

	@param1 data: basicdata object
	@param2 logger: Logging object
	"""
	output = data.geneDir+data.ORFs.split("/")[-1].split(".")[0]+"_prank"
	runPrank(data.ORFs, data.geneName, output)
	data.aln = output+".best.fas"
	logger.info("Performed alignement using Prank")

#######=================================================================================================================
######PhyML=============================================================================================================

def runPhyML(aln, geneDir):
	"""
	Function converting fasta file to phylip and running PhyML.

	@param1 aln: Path
	@param2 geneDir: Gene directory
	@return outPhy: Path to PhyML results file
	"""

	# convert to Phylip format and replace "!" symbols inserted by MACSE for frame shifts by "N"
	outPhy = geneDir+aln.split("/")[-1].split(".")[0]
	tmp = geneDir+aln.split("/")[-1].split(".")[0]+".tmp"
	logger = logging.getLogger("main")
	
	with open(aln, "rU") as aln:
		aln = aln.read().replace("!", "N")
		with open(tmp, "w") as temp:
			temp.write(aln)
			
	input_handle = open(tmp, "rU")
	output_handle = open(outPhy, "w")
	
	alignments = AlignIO.parse(input_handle, "fasta")
	AlignIO.write(alignments, output_handle, "phylip-relaxed")

	output_handle.close()
	input_handle.close()
	os.remove(tmp)

	# PhyML
	cmd("phyml -i {:s} -v e -b -2".format(outPhy), False)
	
	#logger.debug(cmd)
	
	return(outPhy)


def phyMLTree(data, logger):
	"""
	Function creating tree attribute in each gene object of the list.

	@param1 data: List of gene objects
	@param2 logger: Logging object
	"""
	dAlnTree = {}
	logger.info("Beginning of PhyML works")
	TreesFile = runPhyML(data.aln, data.geneDir)
	data.tree = TreesFile+"_phyml_tree.txt"
	logger.info("Reconstructed tree using PhyML")

	dAlnTree[data.aln] = data.tree
	return dAlnTree

#######=================================================================================================================

######GARD==============================================================================================================

def runGARD(aln, o, logger):
	"""
	Function creating the batch file to run GARD (Kosakovsky Pond et al., 2006).
	
	@param1 aln: Path
	@param2 o: Ouput directory
	@return gardRes: Path to GARD output file
	"""

	batchFile = o+aln.split("/")[-1].split(".")[0]+"_gard.bf"
	gardRes = o+aln.split("/")[-1].split(".")[0]+".gard"
	outGard = gardRes+"_out"
	errGard = gardRes+"_err"
	
	with open(batchFile, "w") as bf:
		bf.write("inputRedirect = {};\n")
		bf.write("inputRedirect[\"01\"] = \"{:s}\";\n".format(aln))
		bf.write("inputRedirect[\"02\"] = \"010010\";\n")
		bf.write("inputRedirect[\"03\"] = \"None\";\n")
		bf.write("inputRedirect[\"04\"] = \"{:s}\";\n".format(gardRes))
		bf.write("ExecuteAFile(HYPHY_LIB_DIRECTORY + \"TemplateBatchFiles\" + DIRECTORY_SEPARATOR + \"GARD.bf\", inputRedirect);\n")
	logger.debug("Batch file: {:s}".format(batchFile))
	
	cmd = "mpirun -np 2 HYPHYMPI {:s}".format(batchFile)
	lCmd = shlex.split(cmd)
	with open(outGard, "w") as o, open(errGard, "w") as e:
		runGARD = subprocess.run(lCmd, shell=False, check=True, stdout=o, stderr=e)
	
	logger.debug(cmd)
	return(gardRes)


def procGARD(gardRes, aln):
	"""
	Function creating the batch file to run GARD processor for KH test (Kosakovsky Pond et al., 2006).
	
	@param1 gardRes: Path to GARD output file
	@param2 aln: Path to alignment file
	@return otGardProc: Path to GARDprocessor output file
	"""
	batchFile = gardRes.split(".")[0]+"_gardprocessor.bf"
	splitsFile = gardRes+"_splits"
	outGardProc = gardRes+"_outproc"
	errGardProc = gardRes+"_errproc"
	logger = logging.getLogger("main")
	
	with open(batchFile, "w") as bf:
		bf.write("inputRedirect = {};\n")
		bf.write("inputRedirect[\"01\"] = \"{:s}\";\n".format(aln))
		bf.write("inputRedirect[\"02\"] = \"{:s}\";\n".format(splitsFile))
		bf.write("ExecuteAFile(HYPHY_LIB_DIRECTORY + \"TemplateBatchFiles\" + DIRECTORY_SEPARATOR + \"GARDProcessor.bf\", inputRedirect);\n".format())
	#logger.debug("Batch file: {:s}".format(batchFile))
	
	with open(outGardProc, "w") as o, open(errGardProc, "w") as e:
		runCmd = subprocess.run(["HYPHYMP", batchFile], shell=False, check=True, stdout=o, stderr=e)
	
	return(outGardProc)

def parseGard(kh, aln, geneName, pvalue, o):
	"""
	Function returning the cut fragments following GARD analysis and identification of significant breakpoints.

	@param1 kh: Path to GARDprocessor output file
	@param2 aln: Path to alignment file
	@param3 geneName: Gene name
	@param4 pvalue: Float
	@return lOutFrag: List of Path (Fragments in fasta files)
	"""
	lBP = []
	with open(kh, "r") as f:
		lLine = f.readlines()
	
	index = 0
	while lLine[index].startswith("Breakpoint") != True:
		index += 1

	#If there're breakpoints, add it in lBP
	if lLine[index+1] != "":
		index += 1
		while lLine[index].startswith(" "):
			line = [float(item.strip()) for item in lLine[index].split("|")]
			if line[2] < pvalue and line[4] < pvalue:
				lBP.append(int(line[0]))
			index += 1
	
		print("There are {:d} significant breakpoints in the alignement for gene {:s} at positions {}".format(len(lBP), geneName, lBP))
		
		#If there're breakpoint(s), cut sequence in subsequences according to breakpoints
		if len(lBP) > 0:
			dFname2Fseq = {}
			for fasta in SeqIO.parse(open(aln),'fasta'):
				dFname2Fseq[fasta.id] = str(fasta.seq)
			
			#Creation of a dico where atgc in sequence has been replace by 1 and - by 0
			lSeqBin = []
			lNameGene = []
			for fastaSeq in dFname2Fseq:
				lSeqBin.append(dFname2Fseq[fastaSeq].lower().replace("a", "1").replace("t", "1").replace("c", "1").replace("g", "1").replace("-", "0"))
				lNameGene.append(fastaSeq)

			#looking for a multiple of 3 (number of letter) (subsequence ends on or after the breakpoint)
			nbSeq = len(lNameGene)
			lenSeq = len(lSeqBin[0])
			lPos = [0]
			lBPprec = [0 for i in range(len(lSeqBin))]
			lFrag = []
			for bp in lBP:
				index = 0
				lPosInter = []
				index2 = 0
				for elem in lSeqBin:
					nb = elem[lBPprec[index]:bp].count('1') 
					pos = bp

					while nb%3 != 0:
						pos += 1
						nb = elem[lBPprec[index]:pos].count('1')
					lPosInter.append(pos)
					index2 += 1
				lBPprec = lPosInter
				lPos.append(max(lPosInter))
				lFrag += [ dFname2Fseq[lNameGene[j]][lPos[-2]:lPos[-1]] for j in range(nbSeq) ]

			#Adding subsequences that start at the last breakpoint to the end
			lFrag += [dFname2Fseq[lNameGene[i]][lPos[-1]:] for i in range(nbSeq)]

			lBP = [0]+lBP+[lenSeq]
			lOutFrag = []
			index = 0
			for x in range(1,len(lBP)):
				dFrag = {}
				extension = "{:d}:{:d}".format(lBP[x-1], lBP[x])
				outFrag = o+aln.split("/")[-1].split(".")[0]+"_frag"+extension+".bes.fas"
				for name in lNameGene:
					dFrag[name] = lFrag[index]
					index += 1
				with open(outFrag, "w") as outF:
					outF.write(FastaResFunc.dict2fasta(dFrag))
				lOutFrag.append(outFrag)

			return lOutFrag
		else:
			return []
	else:
		return []

def gardRecomb(data, pvalue, dAT, logger):
	"""
	Procedure which execute gard functions on each alignment given.

	@param1 data: basicData object
	@param2 pvalue: threshold value
	@param3 step: String
	@param4 logger: Logging object
	"""
	try:
		setattr(data, "lGard", [])
		nb = 1	
		for aln in dAT:
			
			directory = data.o+aln.split('/')[-1].split(".")[0]+"/"
			logger.info("trying GARD")
			listGard = runGARD(aln, directory, logger)
				
			logger.info("Checked for recombination using HYPHY GARD")
					
			listGardProc = procGARD(listGard, aln)
			data.lGard.append(listGardProc)
			logger.info("Extracted GARD results")

			listCut = parseGard(listGardProc, aln, data.geneName, float(pvalue), directory)

			logger.info("Parsed GARD results")

			for path in listCut:
				logger.info(path+": Phyml.")
				out = runPhyML(path, directory)+"_phyml_tree.txt"
				dAT[path] = out

			return dAT

	except Exception:
		logger.info("GARD uncountered an unexpected error, skipping.")
		return dAT

#######=================================================================================================================
