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



if __name__ == "__main__" :	

	with open(sys.argv[1], 'r') as config_in:
		config_dict = json.load(config_in)

	parameters = config_dict["parameters"]
	data = config_dict["data"]
	data["firstStep"] = "orf" # a enlever après 
	
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

	config_dict["parameters"] = parameters
	config_dict["data"] = data

	with open(sys.argv[1],'w') as config_out:
		json.dump(config_dict, config_out)


