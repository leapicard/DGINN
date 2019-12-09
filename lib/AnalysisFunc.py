import FastaResFunc, ExtractFunc
import logging, subprocess, shlex, os, ete3
from Bio import SeqIO, AlignIO
from collections import defaultdict, OrderedDict
from itertools import chain
from statistics import median

######Functions=============================================================================================================

def cmd(commandLine, choice):
	"""
	Function executing a command line in a bash terminal.

	@param1 commandLine: String corresponding to a bash command line
	@param2 choice: Boolean determining whether the command is executed within the shell 
	"""
	
	lCmd = shlex.split(commandLine)
	run = subprocess.run(lCmd, 
						 shell=choice, 
						 check=True, 
						 stdout=subprocess.PIPE, 
						 stderr=subprocess.PIPE)

######ORF===================================================================================================================

def getORFs(catFile, queryName, geneDir):
	"""
	Function to find Open Reading Frames within the sequence of each gene and select the longest one.

	@param1 catFile: Path
	@param2 geneName: Gene name
	@param3 geneDir: Gene directory
	@return outORF: Path to the file containing the longest ORFs
	"""

	outORFraw = geneDir+catFile.split("/")[-1].split(".")[0]+"_allORFs.fasta"
	logger = logging.getLogger("main")
	
	cmd("getorf -sequence {:s} -outseq {:s} -table 0 -find 3 -noreverse".format(catFile, outORFraw), False)
	
	logger.debug(cmd)
	
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
				firstOcc = [x for x in dupl if sp in x][0]
			dupl.remove(firstOcc)
		
		for i in dupl:
			dId2Longest.pop(i, None)
			n += 1
			logger.debug("Deleted sequence {:s} (duplicate)".format(i))
		
	logger.info("Deleted {} sequences as duplicates".format(n))
	
	outORF = outORFraw.replace("_allORFs.fasta","_longestORFs.fasta")

	with open(outORF, "w") as outO:
		outO.write(FastaResFunc.dict2fasta(dId2Longest))
	
	logger.info("Extracted longest ORFs: {:s}".format(outORF))

	return(outORF)


def orfFinder(data):
	"""
	Procedure which launch the ORF step

	@param1 data: basicdata object
	@param2 logger: Logging object
	"""
	
	ORFile = getORFs(data.seqFile, 
					 data.queryName, 
					 data.o)
	setattr(data, "ORFs", ORFile)

#######=================================================================================================================
######PRANK=============================================================================================================
	
def runPrank(ORFs, geneName, o):
	"""
	Function to run PRANK software for codon alignment (Löytynoja, 2014).
	
	@param1 ORFs: Path
	@param2 geneName: Gene name
	@param3 geneDir: Gene Directory
	@return outPrank: Path to Prank results file
	"""
	logger = logging.getLogger("main")
	logger.info("Started Prank codon alignment")
	outPrank = o+ORFs.split("/")[-1].split(".")[0]+"_prank"
	
	cmd("prank -d={:s} -o={:s} -codon -F".format(ORFs, outPrank), False)
	
	logger.info("Finished Prank codon alignment: {:s}.best.fas".format(outPrank))
	
	if os.path.exists(outPrank+".fas"):
		os.rename(outPrank+".fas", outPrank+".best.fas")
	
	return(outPrank+".best.fas")

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
	logger = logging.getLogger("main")
	
	if queryName in dId2Seq:
		logger.info("Discarding sequences with less than {:d}% coverage of query.".format(cov))
		outCov = o+aln.split("/")[-1].split(".")[0]+"_mincov.fasta"
		
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
			outC.write(FastaResFunc.dict2fasta(dKeep))
		
		logger.info("Discarded {:d} sequences".format(nbOut))
	
		return(outCov, nbOut)
	
	else:
		logger.warning("Provided query name not found in the alignment, skipping coverage check.")
		return(aln, 0)

def alnPrank(data, logger):
	"""
	Function creating alignment (aln) attribute in each gene object of the list.

	@param1 data: basicdata object
	@param2 logger: Logging object
	"""
	aln = runPrank(data.ORFs, 
				   data.geneName, 
				   data.o)
	data.aln = aln

#######=================================================================================================================
######PhyML=============================================================================================================

def runPhyML(aln, geneDir):
	"""
	Function converting fasta file to phylip and running PhyML.

	@param1 aln: Path
	@param2 geneDir: Gene directory
	@return outPhy: Path to PhyML results file
	"""
	# convert to Phylip format and replace eventual "!" symbols (relic from using MACSE)
	outPhy = geneDir+aln.split("/")[-1].split(".")[0]+".phylip"
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
	logger.debug("phyml -i {:s} -v e -b -2".format(outPhy))
	
	return(outPhy)


def phyMLTree(data, logger):
	"""
	Function creating tree attribute in each gene object of the list.

	@param1 data: List of gene objects
	@param2 logger: Logging object
	@return dAlnTree: Updated dictionary of alignments and their corresponding trees
	"""
	dAlnTree = {}
	logger.info("Running PhyML to produce gene phylogenetic tree")
	TreesFile = runPhyML(data.aln, data.o)
	data.tree = TreesFile+"_phyml_tree.txt"
	logger.info("Reconstructed tree using PhyML: {:s}".format(data.tree))

	dAlnTree[data.aln] = data.tree
	return(dAlnTree)

def cutLongBranches(aln, dAlnTree, logger):
	"""
	Check for overly long branches in a tree and separate both tree and corresponding alignment if found.
	
	@param1 aln: Fasta alignment
	@param2 tree: Tree corresponding to the alignment
	@param3 logger: Logging object
	@return dAlnTree: Updated dictionary of alignments and their corresponding trees
	"""
	logger.info("Looking for long branches.")
	loadTree = ete3.Tree(dAlnTree[aln])
	dist = [leaf.dist for leaf in loadTree.traverse()]
	medianDist = median(dist)
	longDist = medianDist * 50
	matches = [leaf for leaf in loadTree.traverse() if leaf.dist>50.0]
	
	if len(matches) > 0:
		logger.info("{} long branches found, separating alignments.".format(len(matches)))
		
		seqs = SeqIO.parse(open(aln),'fasta')
		dID2Seq = {gene.id: gene.seq for gene in seqs}
		
		for node in matches:
			gp = node.get_children()
			lNewGp = list(chain.from_iterable([x.get_leaf_names() for x in gp]))
		
			newAln = aln.split(".")[0]+"_split"+str(matches.index(node)+1)+".fasta"
			
			dNewAln = {gene:dID2Seq[gene] for gene in lNewGp}
			for k in lNewGp:
				dID2Seq.pop(k, None)
			
			# create new file of sequences
			with open(newAln, "w") as fasta:
				fasta.write(FastaResFunc.dict2fasta(dNewAln))
			
			dAlnTree[newAln] = ""
		
		alnLeft = aln.split(".")[0]+"_split"+str(len(matches)+1)+".fasta"
		with open(alnLeft, "w") as fasta:
			fasta.write(FastaResFunc.dict2fasta(dID2Seq))
		
		dAlnTree[alnLeft] = ""
		dAlnTree.pop(aln, None)
		
	else:
		logger.info("No long branches found.")
	
	return(dAlnTree)

def checkPhyMLTree(data, dAlnTree, logger):
	dAlnTree = cutLongBranches(data.aln, dAlnTree, logger)
	dAlnTree2 = {}
	
	for aln in dAlnTree:
		if dAlnTree[aln] == "":
			logger.info("Reconstructing alignments and phylogenies following long branch parsing.")
			aln = runPrank(aln, data.geneName, data.o)
			tree = runPhyML(aln, data.o)
			dAlnTree2[aln] = tree+"_phyml_tree.txt"
			
	dAlnTree.update(dAlnTree2)
		
	return(dAlnTree)
		
#######=================================================================================================================

######GARD==============================================================================================================

def runGARD(aln, o, hostFile, logger):
	"""
	Function creating the batch file to run GARD (Kosakovsky Pond et al., 2006).
	
	@param1 aln: Path
	@param2 o: Ouput directory
	@return gardRes: Path to GARD output file
	"""

	gardRes = o+aln.split("/")[-1].split(".")[0]+".gard"
	gardJson = o+aln.split("/")[-1].split(".")[0]+".json"
	outGard = gardRes+"_out"
	errGard = gardRes+"_err"
	
	# valid for hyphy 2.3
	batchFile = o+aln.split("/")[-1].split(".")[0]+"_gard.bf"
	with open(batchFile, "w") as bf:
		bf.write("inputRedirect = {};\n")
		bf.write("inputRedirect[\"01\"] = \"{:s}\";\n".format(aln))
		bf.write("inputRedirect[\"02\"] = \"010010\";\n")
		bf.write("inputRedirect[\"03\"] = \"None\";\n")
		bf.write("inputRedirect[\"04\"] = \"{:s}\";\n".format(gardRes))
		bf.write("ExecuteAFile(HYPHY_LIB_DIRECTORY + \"TemplateBatchFiles\" + DIRECTORY_SEPARATOR + \"GARD.bf\", inputRedirect);\n")
	logger.debug("Batch file: {:s}".format(batchFile))
	
	"""
	#for hyphy 2.5
	if hostFile == "":
		cmd = "mpirun -np 2 HYPHYMPI GARD --type codon --model HKY --alignment {:s} --output {:s} --output-lf {:s}".format(aln,
																														   gardRes,
																														   gardJson)
	else:
		cmd = "mpirun -np 2 -hostfile {:s} HYPHYMPI GARD --type codon --model HKY --alignment {:s} --output {:s} --output-lf {:s}".format(hostfile,
																																		  aln,
																									  									  gardRes,
																									  									  gardJson)
	"""

	if hostFile == "":
		cmd = "mpirun -np 2 HYPHYMPI {:s}".format(batchFile)
	else:
		cmd = "mpirun -np 2 -hostfile {:s} HYPHYMPI {:s}".format(hostFile, batchFile)


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
	#logger = logging.getLogger("main")
	
	with open(batchFile, "w") as bf:
		bf.write("inputRedirect = {};\n")
		bf.write("inputRedirect[\"01\"] = \"{:s}\";\n".format(aln))
		bf.write("inputRedirect[\"02\"] = \"{:s}\";\n".format(splitsFile))
		bf.write("ExecuteAFile(HYPHY_LIB_DIRECTORY + \"TemplateBatchFiles\" + DIRECTORY_SEPARATOR + \"GARDProcessor.bf\", inputRedirect);\n".format())
	#logger.debug("Batch file: {:s}".format(batchFile))
	
	cmd = "HYPHYMP {:s}".format(batchFile)
	lCmd = shlex.split(cmd)
	with open(outGardProc, "w") as o, open(errGardProc, "w") as e:
		runGARD = subprocess.run(lCmd, shell=False, check=True, stdout=o, stderr=e)
	
	return(outGardProc)

def parseGard(kh, aln, pvalue, o, logger):
	"""
	Function returning the cut fragments following GARD analysis and identification of significant breakpoints.

	@param1 kh: Path to GARDprocessor output file
	@param2 aln: Path to alignment file
	@param3 pvalue: Float
	@param4 o: Path to output directory
	@return lOutFrag: List of Path (Fragments in fasta files)
	"""
	lBP = []
	with open(kh, "r") as f:
		lLine = f.readlines()
		finalIndex = len(lLine)
	
	index = 0
	while lLine[index].startswith("Breakpoint") != True and index < finalIndex:
		index += 1
	
	#If there are breakpoints, add it in lBP
	if lLine[index+1] != "":
		index += 1
		while lLine[index].startswith(" "):
			line = [float(item.strip()) for item in lLine[index].split("|")]
			if line[2] < pvalue and line[4] < pvalue:
				lBP.append(int(line[0]))
			index += 1
		
		if len(lBP) > 0:
			logger.info("There are {:d} significant breakpoints in alignement {:s} at positions {}".format(len(lBP), aln, lBP))
		else:
			logger.info("There are no significant breakpoints in alignement {:s}.".format(aln))
		
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
				while bp%3 != 0:
					bp += 1
				lPos.append(bp)
				lFrag += [ dFname2Fseq[lNameGene[j]][lPos[-2]:lPos[-1]] for j in range(nbSeq) ]
			
			#Adding subsequences that start at the last breakpoint to the end
			lFrag += [dFname2Fseq[lNameGene[i]][lPos[-1]:] for i in range(nbSeq)]

			lBP = lPos+[lenSeq]
			lOutFrag = []
			index = 0
			for x in range(1,len(lBP)):
				dFrag = {}
				if lBP[x-1] == 0:
					extension = "{:d}to{:d}".format(lBP[x-1], lBP[x])
				else:
					extension = "{:d}to{:d}".format(lBP[x-1]-1, lBP[x])

				outFrag = o+aln.split("/")[-1].split(".")[0]+"_frag"+extension+".best.fas"
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

def gardRecomb(data, pvalue, dAT, hostFile, logger):
	"""
	Procedure which execute gard functions on each alignment given.

	@param1 data: basicData object
	@param2 pvalue: threshold value
	@param3 step: String
	@param4 logger: Logging object
	"""
	#try:
	setattr(data, "lGard", [])
	nb = 1
	dFrag = {}
	for aln in dAT:
		logger.info("Running GARD on {:s}".format(aln))
		gardRes = runGARD(aln, 
						  data.o, 
						  hostFile, 
						  data.logger)
			
		logger.info("Checked for recombination using HYPHY GARD.")
				
		procGardRes = procGARD(gardRes, aln)
		#data.lGard.append(listGardProc)
		logger.info("Extracted GARD results.")

		parsedFragments = parseGard(procGardRes, 
									aln, 
									float(pvalue), 
									data.o, 
									logger)

		logger.info("Parsed GARD results.")
		if len(parsedFragments) > 0:
			for frag in parsedFragments:
				logger.info("Running Phyml on fragment: {:s}".format(frag))
				fragTree = runPhyML(frag, data.o)+"_phyml_tree.txt"
				dFrag[frag] = fragTree
		else:
			logger.info("No fragments of recombination identified.")
			
	dAT.update(dFrag)
	return(dAT)
	
	"""except Exception:
		logger.info("GARD uncountered an unexpected error, skipping.")
		return dAT"""

#######=================================================================================================================
