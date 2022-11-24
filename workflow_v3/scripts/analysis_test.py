import sys
import logging, subprocess, shlex, os, ete3
from Bio import SeqIO, AlignIO
from collections import defaultdict, OrderedDict
from itertools import chain
from statistics import median, mean
import pandas as pd
import json
import fastares_test, loadfile_test, tree_test


######Functions=============================================================================================================

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


def runPhyML(aln, phymlOpt, geneDir):
	"""
	Function converting fasta file to phylip and running PhyML.

	@param1 aln: Path
	@param2 geneDir: Gene directory
	@return outPhy: Path to PhyML results file
	"""
	# convert to Phylip format and replace eventual "!" symbols (relic from using MACSE)
	origin = os.getcwd()
	os.chdir(geneDir)
	outPhy = aln.split("/")[-1].split(".")[0]+".phylip"
	tmp = aln.split("/")[-1].split(".")[0]+".tmp"
	aln = aln.split("/")[-1]

	logger = logging.getLogger("main.tree")
	
	with open(aln, "r") as aln2:
		laln = aln2.read().replace("!", "N")
		aln2.close()
		with open(tmp, "w") as temp:
			temp.write(laln)
			temp.close()

	input_handle = open(tmp, "r")
	output_handle = open(outPhy, "w")
	
	alignments = AlignIO.parse(input_handle, "fasta")
	AlignIO.write(alignments, output_handle, "phylip-relaxed")

	output_handle.close()
	input_handle.close()
	os.remove(tmp)
	
	# PhyML
	if phymlOpt != "":
		try:
			opt=phymlOpt.split("ALN ")[1]
			logger.debug("phyml -i {:s} {}".format(outPhy, opt))
			cmd("phyml -i {:s} {}".format(outPhy, opt), False)
		except:
			logger.info("PhyML couldn't run with the provided info {}, running with default options.".format(phymlOpt))
			cmd("phyml -i {:s} -v e -b -2".format(outPhy), False)
	else:
		logger.debug("phyml -i {:s} -v e -b -2".format(outPhy))
		cmd("phyml -i {:s} -v e -b -2".format(outPhy), False)
		
	os.chdir(origin)
	
	return(geneDir+outPhy)


def phyMLTree(data,phymlOpt):
	"""
	Function creating tree attribute in each gene object of the list.

	@param1 data: List of gene objects
	@return dAlnTree: Updated dictionary of alignments and their corresponding trees
	"""

	logger = logging.getLogger("main.tree")
	dAlnTree = {}
	logger.info("Running PhyML to produce gene phylogenetic tree")
	TreesFile = runPhyML(data["aln"], phymlOpt, data["o"])
	data["tree"] = TreesFile+"_phyml_tree.txt"
	logger.info("Reconstructed tree using PhyML: {:s}".format(data["tree"]))
	dAlnTree[data["aln"]] = data["tree"]
	return dAlnTree
    
def cutLongBranches(aln, dAlnTree, nbSp, LBOpt, logger):
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
	#longDist = 500
	
	if "cutoff" in LBOpt:
		if "(" in LBOpt:
			factor = float(LBOpt.split("(")[1].replace(")", ""))
		else:
			factor = 50
		medianDist = median(dist)
		meanDist = mean(dist)
		longDist = meanDist * factor
	elif "IQR" in LBOpt:
		if "(" in LBOpt:
			factor = float(LBOpt.split("(")[1].replace(")", ""))
		else:
			factor = 50
		df = pd.DataFrame(dist)
		Q1 = df.quantile(0.25)
		Q3 = df.quantile(0.75)
		IQR = Q3 - Q1
		lDist = Q3 + (factor * IQR)
		longDist = lDist[0]
	
	logger.info("Long branches will be evaluated through the {} method (factor {})".format(LBOpt, factor))
	nbSp = int(nbSp)	
	matches = [leaf for leaf in loadTree.traverse() if leaf.dist>longDist]

	if len(matches) > 0:
		logger.info("{} long branches found, separating alignments.".format(len(matches)))
		
		seqs = SeqIO.parse(open(aln),'fasta')
		dID2Seq = {gene.id: gene.seq for gene in seqs}
				
		for node in matches:
			gp = node.get_children()
			lNewGp = list(chain.from_iterable([x.get_leaf_names() for x in gp]))
		
			newAln = aln.split(".")[0]+"_part"+str(matches.index(node)+1)+".fasta"
			
			dNewAln = {gene:dID2Seq[gene] for gene in lNewGp if gene in dID2Seq}
			for k in lNewGp:
				dID2Seq.pop(k, None)
			
			# create new file of sequences
			
			if len(dNewAln) > nbSp - 1:
				with open(newAln, "w") as fasta:
					fasta.write(fastares_test.dict2fasta(dNewAln))
					fasta.close()		  
				dAlnTree[newAln] = ""
			else:
				logger.info("Sequences {} will not be considered for downstream analyses as they do not compose a large enough group.".format(dNewAln.keys()))
			
		alnLeft = aln.split(".")[0]+"_part"+str(len(matches)+1)+".fasta"

		if len(dID2Seq) > nbSp - 1:
			with open(alnLeft, "w") as fasta:
				fasta.write(fastares_test.dict2fasta(dID2Seq))
				logger.info("\tNew alignment:%s"%{alnLeft})
				fasta.close()	  
			dAlnTree[alnLeft] = ""
		else:
			logger.info("Sequences in {} will not be considered for downstream analyses as they do not compose a large enough group.".format(dID2Seq.keys()))
			
		dAlnTree.pop(aln, None)
		
	else:
		logger.info("No long branches found.")
		
	return(dAlnTree)

def checkPhyMLTree(data, dAlnTree, nbSp, LBopt, step="duplication"):
	logger=logging.getLogger(".".join(["main",step]))
	dAlnTree = cutLongBranches(data["aln"], dAlnTree, nbSp, LBopt, logger)
	dAlnTree2 = {}
	
	for aln in dAlnTree:
		if dAlnTree[aln] == "":
			logger.info("Reconstructing alignments and phylogenies following long branch parsing.")
			#aln = runPrank(aln, data["o"]) #pb à régler 
			tree = runPhyML(aln, LBopt, data["o"])
			dAlnTree2[aln] = tree+"_phyml_tree.txt"
			
	#dAlnTree.update(dAlnTree2)
	if len(dAlnTree2) > 0:
		return(dAlnTree2)
	else:
		return(dAlnTree)

def runGARD(aln, o, hostFile, logger):
	"""
	Function creating the batch file to run GARD (Kosakovsky Pond et al., 2006).
	
	@param1 aln: Path
	@param2 o: Ouput directory
	@return gardRes: Path to GARD output file
	"""

	gardRes = o+aln.split("/")[-1].split(".")[0]+".gard"
	gardJson = o+aln.split("/")[-1]+".GARD.json"
	outGard = gardRes+"_out"
	errGard = gardRes+"_err"
	
	#for hyphy 2.5
	if hostFile == "":
		cmd = "mpirun -np 4 HYPHYMPI GARD --alignment {:s} --output {:s} --output-lf {:s}".format(aln, gardJson, gardRes)
	else:
		cmd = "mpirun -np 4 -hostfile {:s} HYPHYMPI GARD --alignment {:s} --output {:s} --output-lf {:s}".format(hostFile, aln, gardJson, gardRes)
          
	lCmd = shlex.split(cmd)
	with open(outGard, "w") as o, open(errGard, "w") as e:
		runGARD = subprocess.run(lCmd, shell=False, check=True, stdout=o, stderr=e)
		o.close()
		e.close()
	logger.debug(cmd)
	logger.info(gardJson)
	return(gardJson)



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
	logger = logging.getLogger("main.recombination")
	
	with open(batchFile, "w") as bf:
		bf.write("inputRedirect = {};\n")
		bf.write("inputRedirect[\"01\"] = \"{:s}\";\n".format(aln))
		bf.write("inputRedirect[\"02\"] = \"{:s}\";\n".format(splitsFile))
		bf.write("ExecuteAFile(HYPHY_LIB_DIRECTORY + \"TemplateBatchFiles\" + DIRECTORY_SEPARATOR + \"GARDProcessor.bf\", inputRedirect);\n".format())
	#logger.debug("Batch file: {:s}".format(batchFile))
	
	cmd = "hyphy {:s}".format(batchFile)
	logger.info(cmd)
	logger.info(os.system(cmd))
	logger.info("====================")
	lCmd = shlex.split(cmd)
	with open(outGardProc, "w") as o, open(errGardProc, "w") as e:
		runGARD = subprocess.run(lCmd, shell=False, check=True, stdout=o, stderr=e)
		o.close()
		e.close()
	return(outGardProc)

def parseGard(kh, aln, o, logger):
	"""
	Function returning the cut fragments following GARD analysis and identification of significant breakpoints.

	@param1 kh: Path to GARD.json output file
	@param2 aln: Path to alignment file
	@param3 pvalue: Float
	@param4 o: Path to output directory
	@return lOutFrag: List of Path (Fragments in fasta files)
	"""
	lBP = []
	f = open(kh, "r")
	lLine = f.readline()
	while lLine:
		if lLine.find("\"breakpoints\"")!=-1:
			lLine = f.readline()
			lLine=lLine[lLine.find("[")+1:lLine.find("]")]
			lBP=list(map(int, lLine.split(",")))
			break
		lLine = f.readline()
	f.close()
	index = 0
	
	#If there are breakpoints, add it in lBP
	if len(lBP) > 0:
		logger.info("There are {:d} significant breakpoints in alignement {:s} at positions {}".format(len(lBP), aln, lBP))
	else:
		logger.info("There are no significant breakpoints in alignement {:s}.".format(aln))
		return []
              
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
				extension = "_{:d}_{:d}".format(lBP[x-1], lBP[x])
			else:
				extension = "_{:d}_{:d}".format(lBP[x-1]-1, lBP[x])

			outFrag = o+aln.split("/")[-1].split(".")[0]+"_frag"+extension+".best.fas"
			for name in lNameGene:
				dFrag[name] = lFrag[index]
				index += 1
			with open(outFrag, "w") as outF:
				outF.write(fastares_test.dict2fasta(dFrag))
				logger.info("\tNew alignment: %s"%{outFrag})
				outF.close()
				lOutFrag.append(outFrag)

		return lOutFrag
	else:
		return []

def gardRecomb(data, dAT, hostFile):
	"""
	Procedure which execute gard functions on each alignment given.

	@param1 data: basicData object
	@param2 pvalue: threshold value
	@param3 step: String
	@param4 logger: Logging object
	"""
	#try:
	logger=logging.getLogger("main.recombination")
	data["lGard"] = []
	nb = 1
	dFrag = {}
	for aln in dAT:
		logger.info("Running GARD on {:s}".format(aln))
		gardRes = runGARD(aln, 
			    data["o"], 
			    hostFile, 
			    data["logger"])
			
		logger.info("Checked for recombination using HYPHY GARD.")
				
			
		parsedFragments = parseGard(gardRes, aln, data["o"], logger)

		logger.info("Parsed GARD results.")
		if len(parsedFragments) > 0:
			for frag in parsedFragments:
				logger.info("Running Phyml on fragment: {:s}".format(frag))
				fragTree = runPhyML(frag, "",  data["o"])
				dFrag[frag] = fragTree+"_phyml_tree.txt"
		else:
			logger.info("No fragments of recombination identified.")
	      
		dAT = dFrag
		return(dAT)
	


if __name__ == "__main__" :	

	with open(sys.argv[1], 'r') as config_in:
		config_dict = json.load(config_in)

	parameters = config_dict["parameters"]
	data = config_dict["data"]
	data["firstStep"] = "orf" # a enlever après 
	parameters["hostfile"] = "" # a enlever après 
	data["logger"] = logging.getLogger("main")

	if sys.argv[2] == "phyMLTree":

		if parameters["step"] == "tree":
			data = loadfile_test.phymlRecEntry(data)

			if "tree" == data["firstStep"]:
				loadfile_test.spTreeCheck(data, 
							data["firstStep"], 
							parameters["duplication"])

		dAlTree = phyMLTree(data, parameters["phymlOpt"])
		data["dAlTree"] = dAlTree


	elif sys.argv[2] == "checkPhyMLTree":
		if parameters["duplication"]:
			if parameters["step"] == "duplication":
				data, data["dAlTree"] = loadfile_test.duplPSEntry(data)

				if data["firstStep"] == "duplication":
					dTree = data["dAlTree"].pop(data["aln"])
					loadfile_test.spTreeCheck(data, 
								data["firstStep"], 
								parameters["duplication"])
					data["dAlTree"][data["aln"]] = dTree
							
			dAlTree = checkPhyMLTree(data, data["dAlTree"], parameters["nbspecies"], parameters["LBopt"])
			dAlTree = tree_test.treeTreatment(data, data["dAlTree"], parameters["nbspecies"], parameters["phymlOpt"])
			data["dAlTree"] = dAlTree


	elif sys.argv[2] == "gardRecomb":
		if parameters["recombination"] :
			if parameters["step"] == "recombination":
				data = loadfile_test.phymlRecEntry(data, "main.recombination")
				data["dAlTree"][data["aln"]] = ""

		data["dAlTree"] = gardRecomb(data, data["dAlTree"], parameters["hostfile"])


	config_dict["parameters"] = parameters
	config_dict["data"] = data

	with open(sys.argv[1],'w') as config_out:
		json.dump(config_dict, config_out)


