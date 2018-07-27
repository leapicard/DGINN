from scipy import stats
from collections import OrderedDict
import re, os, logging

def getParams(models, bppml, mixed, opbFile, gnhFile, Busted, Meme, opb, gnh):
	# Check analyses to be run and where the parameters file are

	dCtrls = {}
	lModels = []
	if models != "":
		dCtrls["bppml"] = bppml
		lModels = re.compile("\s*,\s*").split(models)
		dCtrls["bppmixedlikelihood"] = mixed
	if opbFile != "" and opb == "True":
		dCtrls["OPB"] = opbFile
	if gnh != "" and gnh == "True":
		dCtrls["GNH"] = gnhFile
	if Busted == "True":
		dCtrls["BUSTED"] = "/usr/local/lib/hyphy/TemplateBatchFiles/SelectionAnalyses/BUSTED.bf"
	if Meme == "True":
		dCtrls["MEME"] = "/usr/local/lib/hyphy/TemplateBatchFiles/SelectionAnalyses/MEME.bf"

	return dCtrls, lModels

def supBoot(outDir, baseName, treeFile, logger):
	# Suppress bootstrap numbers from treeFile (necessary for HYPHY)
	cladoFile = outDir+baseName+"_clado.tree"
	with open(cladoFile, "w") as clado:
		clado.write(leanTree(treeFile))
	
	logger.debug(leanTree(treeFile))
	return cladoFile

def nbNode(treeFile, outDir, logger):
	# count number of nodes in tree file
	with open(treeFile, "r") as tree:
		data = tree.read()
		nodes = str(data.count("(")+data.count(")"))
		logger.info("There are %s nodes in the provided tree." %nodes)
	logger.info("Output folder: %s" %outDir)

	return nodes

def LRT(ll1, ll2, df):
	"""
	Calculates likelihood ratio test between two models.
	:params ll1, ll2: likelihood of the two models studied
	:param df: degrees of freedom of difference between the two models
	"""
	LR = abs(2*(ll1-ll2))
	stats.chisqprob = lambda chisq, df: stats.chi2.sf(LR, df)
	p = stats.chisqprob(LR, df)
	return(LR, p)

def NHXTree(tree):
	"""
	Take a newick tree and outputs it in NHX format with branch numbering.
	"""
	pattern = re.compile(r"[,\)]")
	par=pattern.findall(tree)
	lm=pattern.split(tree)
	
	nhxpar=["[&&NHX:ND=%d]"%i+par[i] for i in range(len(par))]
	
	lc=["".join(x) for x in zip(lm,nhxpar)]
	sout="".join(lc)+";"
	
	return sout
	
def leanTree(tree):
	"""
	Take a newick tree and returns the cladogram
	"""
	
	pattern = re.compile(r":\d+\.\d+|\n")
	pattern2 = re.compile(r"\)\d+\.\d+")
	#pattern3 = re.compile("\|.+\.\d")
	with open(tree, "r") as tf:
		tf = re.sub(pattern, "", tf.read())
		tf = re.sub(pattern2, ")", tf)
		tf = re.sub("\||\.", "", tf)
		return(tf)
	
def pspFileCreation(path, option):
	if option == "bppml":
		with open(path, "w") as bppml:
			bppml.write("alphabet = Codon(letter=DNA)\ninput.sequence.file=$(INPUTFILE)\ninput.sequence.format=$(FORMAT)\ninput.sequence.sites_to_use = all\ninput.sequence.max_gap_allowed = 50%\ninput.sequence.max_unresolved_allowed = 100%\ninput.tree.file = $(TREEFILE)\ninput.tree.format = Newick\nnonhomogeneous = general\nnonhomogeneous.number_of_models = 1\nnonhomogeneous.root_freq=F3X4(initFreqs=observed)\nmodel1 = $(MODEL)\nmodel1.nodes_id=0:$(NODES)\nlikelihood.recursion = simple\nlikelihood.recursion_simple.compression = recursive\noptimization = FullD(derivatives=Newton)\noptimization.ignore_parameters = $(IGNORE)\noptimization.max_number_f_eval = 100\noptimization.tolerance = 0.01\noutput.tree.file = $(OUTTREE)\noutput.tree.format = Newick\noutput.estimates = $(OUTPARAMS)\noptimization.backup.file = $(BACKUP)")
	if option == "mixedlikelihood":
		with open(path, "w") as mixed:
			mixed.write("alphabet = Codon(letter=DNA)\ninput.sequence.file=$(INPUTFILE)\ninput.sequence.format=$(FORMAT)\ninput.sequence.sites_to_use = all\ninput.sequence.max_gap_allowed = 50%\ninput.sequence.max_unresolved_allowed = 100%\ninput.tree.file = $(TREEFILE)\ninput.tree.format = Newick\nparams = $(PARAMS)\noutput.likelihoods.file = $(OUTINFO)")
	if option == "opbFile":
		with open(path, "w") as opbFile:
			opbFile.write("alphabet = Codon(letter=DNA)\ninput.sequence.file=$(INPUTFILE)\ninput.sequence.format=$(FORMAT)\ninput.sequence.sites_to_use = all\ninput.sequence.max_gap_allowed = 50%\ninput.sequence.max_unresolved_allowed = 100%\ninput.tree.file = $(TREEFILE)\ninput.tree.format = Newick\nnonhomogeneous = one_per_branch\nnonhomogeneous.root_freq=F3X4(initFreqs=observed)\nmodel = YNGP_M0(frequencies=F3X4(initFreqs=observed))\nnonhomogeneous_one_per_branch.shared_parameters = YN98.kappa, YN98.*theta*\nlikelihood.recursion = simple\nlikelihood.recursion_simple.compression = recursive\noptimization = FullD(derivatives=Newton)\noptimization.ignore_parameters = $(IGNORE)\noptimization.max_number_f_eval = 10\noptimization.tolerance = 0.01\noutput.tree.file = $(OUTTREE)\noutput.tree.format = Newick\noutput.estimates = $(OUTPARAMS)\noptimization.backup.file = $(BACKUP)")
	if option == "gnhFile":
		with open(path, "w") as gnhFile:
			gnhFile.write("alphabet = Codon(letter=DNA)\ninput.sequence.file=$(INPUTFILE)\ninput.sequence.format=$(FORMAT)\ninput.sequence.sites_to_use = all\ninput.sequence.max_gap_allowed = 50%\ninput.sequence.max_unresolved_allowed = 100%\ninput.tree.file = $(TREEFILE)\ninput.tree.format = Newick\nnonhomogeneous = general\nnonhomogeneous.number_of_models = 2\nnonhomogeneous.root_freq=F3X4(initFreqs=observed)\nmodel1 = YNGP_M1(frequencies=F3X4(initFreqs=observed))\nmodel1.nodes_id=$(NODE)\nmodel2 = YNGP_M2(frequencies=F3X4(initFreqs=observed),kappa=YNGP_M1.kappa_1,omega=YNGP_M1.omega_1)\nmodel2.nodes_id=$(NODE)\nlikelihood.recursion = simple\nlikelihood.recursion_simple.compression = recursive\noptimization = FullD(derivatives=Newton)\noptimization.ignore_parameters = BrLen,YNGP_M0*,*kappa*,*theta*,Ancient\noptimization.max_number_f_eval = 1000\noptimization.tolerance = 0.0001\noutput.tree.file = $(OUTTREE)\noutput.tree.format = Newick\noutput.estimates = $(OUTPARAMS)\noptimization.backup.file = $(BACKUP)")